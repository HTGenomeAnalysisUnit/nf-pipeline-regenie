#!/usr/bin/env bash
# =============================================================================
# run-tests-wdl.sh
# =============================================================================
# Run WDL test configurations corresponding to the Nextflow run-tests.sh
#
# Usage:
#   ./run-tests-wdl.sh [engine] [backend] [test_name]
#
# Arguments:
#   engine   - cromwell | miniwdl               (default: cromwell)
#   backend  - local | slurm | lsf              (default: local)
#   test_name - optional specific test to run    (e.g., test1)
#
# Examples:
#   ./run-tests-wdl.sh cromwell local           # All tests, Cromwell, local
#   ./run-tests-wdl.sh cromwell slurm           # All tests, Cromwell, SLURM
#   ./run-tests-wdl.sh cromwell lsf test1       # test1 only, Cromwell, LSF
#   ./run-tests-wdl.sh miniwdl local            # All tests, miniwdl, local
#   ./run-tests-wdl.sh miniwdl slurm test7      # test7, miniwdl, SLURM
#
# Environment variables:
#   CROMWELL_JAR      - Path to cromwell.jar (default: cromwell.jar)
#   SLURM_PARTITION   - SLURM partition (optional)
#   SLURM_ACCOUNT     - SLURM account (optional)
#   LSF_QUEUE         - LSF queue name (optional)
#   LSF_PROJECT       - LSF project name (optional)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WDL_DIR="${SCRIPT_DIR}"
TESTS_DIR="${WDL_DIR}/tests"
MAIN_WDL="${WDL_DIR}/regenie_gwas.wdl"
OPTIONS_JSON="${WDL_DIR}/options.json"
BACKENDS_DIR="${WDL_DIR}/backends"
OUTPUT_DIR="${SCRIPT_DIR}/../tests/output_wdl"

# Default engine and backend
ENGINE="${1:-cromwell}"
BACKEND="${2:-local}"
# Optional: run only a specific test
SPECIFIC_TEST="${3:-}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Resolve backend config file for Cromwell
BACKEND_CONF=""
case "$BACKEND" in
    local)
        BACKEND_CONF="${BACKENDS_DIR}/local.conf"
        ;;
    slurm)
        BACKEND_CONF="${BACKENDS_DIR}/slurm.conf"
        ;;
    lsf)
        BACKEND_CONF="${BACKENDS_DIR}/lsf.conf"
        ;;
    *)
        echo -e "${RED}Unknown backend: ${BACKEND}. Use 'local', 'slurm', or 'lsf'.${NC}"
        exit 1
        ;;
esac

echo "=============================================="
echo " REGENIE WDL Pipeline - Test Runner"
echo "=============================================="
echo "Engine:    ${ENGINE}"
echo "Backend:   ${BACKEND}"
echo "WDL file:  ${MAIN_WDL}"
echo "Tests dir: ${TESTS_DIR}"
echo "Output:    ${OUTPUT_DIR}"
echo "=============================================="

# Resolve all paths in JSON to absolute paths
# (WDL engines typically need absolute paths, not relative)
resolve_paths_in_json() {
    local json_file="$1"
    local base_dir
    base_dir="$(cd "$(dirname "$json_file")" && pwd)"
    
    # Replace ../ relative paths with absolute paths
    sed "s|\"\\.\\.\/|\"${base_dir}/../|g" "$json_file" | \
    python3 -c "
import json, sys, os
data = json.load(sys.stdin)
base = '${base_dir}'
resolved = {}
for k, v in data.items():
    if k.startswith('##'):
        continue
    if isinstance(v, str) and v.startswith(base):
        resolved[k] = os.path.realpath(v)
    elif isinstance(v, list):
        new_list = []
        for item in v:
            if isinstance(item, dict):
                new_item = {}
                for ik, iv in item.items():
                    if isinstance(iv, str) and iv.startswith(base):
                        new_item[ik] = os.path.realpath(iv)
                    else:
                        new_item[ik] = iv
                new_list.append(new_item)
            elif isinstance(item, str) and item.startswith(base):
                new_list.append(os.path.realpath(item))
            else:
                new_list.append(item)
        resolved[k] = new_list
    else:
        resolved[k] = v
json.dump(resolved, sys.stdout, indent=4)
"
}

run_test_cromwell() {
    local test_name="$1"
    local inputs_json="$2"
    local resolved_json="${OUTPUT_DIR}/${test_name}_resolved_inputs.json"
    
    echo -e "${YELLOW}>>> Running ${test_name} with Cromwell (${BACKEND})...${NC}"
    
    # Resolve relative paths to absolute
    resolve_paths_in_json "$inputs_json" > "$resolved_json"
    
    java -Dconfig.file="${BACKEND_CONF}" \
        -jar "${CROMWELL_JAR:-cromwell.jar}" run \
        "${MAIN_WDL}" \
        --inputs "$resolved_json" \
        --options "${OPTIONS_JSON}" \
        2>&1 | tee "${OUTPUT_DIR}/${test_name}_cromwell.log"
    
    local exit_code=${PIPESTATUS[0]}
    if [[ $exit_code -eq 0 ]]; then
        echo -e "${GREEN}>>> ${test_name}: PASSED${NC}"
    else
        echo -e "${RED}>>> ${test_name}: FAILED (exit code: ${exit_code})${NC}"
    fi
    return $exit_code
}

run_test_miniwdl() {
    local test_name="$1"
    local inputs_json="$2"
    local resolved_json="${OUTPUT_DIR}/${test_name}_resolved_inputs.json"
    
    echo -e "${YELLOW}>>> Running ${test_name} with miniwdl (${BACKEND})...${NC}"
    
    # Resolve relative paths to absolute
    resolve_paths_in_json "$inputs_json" > "$resolved_json"

    # Set miniwdl config for SLURM if needed
    local miniwdl_env=()
    if [[ "$BACKEND" == "slurm" ]]; then
        miniwdl_env=(env MINIWDL__CFG="${BACKENDS_DIR}/miniwdl.cfg")
    elif [[ "$BACKEND" == "lsf" ]]; then
        echo -e "${RED}NOTE: miniwdl does not natively support LSF. Use Cromwell for LSF.${NC}"
        echo -e "${YELLOW}Falling back to local execution...${NC}"
        miniwdl_env=(env)
    fi
    
    "${miniwdl_env[@]}" miniwdl run \
        "${MAIN_WDL}" \
        --input "$resolved_json" \
        --dir "${OUTPUT_DIR}/${test_name}" \
        2>&1 | tee "${OUTPUT_DIR}/${test_name}_miniwdl.log"
    
    local exit_code=${PIPESTATUS[0]}
    if [[ $exit_code -eq 0 ]]; then
        echo -e "${GREEN}>>> ${test_name}: PASSED${NC}"
    else
        echo -e "${RED}>>> ${test_name}: FAILED (exit code: ${exit_code})${NC}"
    fi
    return $exit_code
}

# Track results
PASSED=0
FAILED=0
SKIPPED=0

for inputs_json in "${TESTS_DIR}"/test*_inputs.json; do
    test_name=$(basename "$inputs_json" _inputs.json)
    
    # If a specific test was requested, skip others
    if [[ -n "$SPECIFIC_TEST" && "$test_name" != "$SPECIFIC_TEST" ]]; then
        continue
    fi
    
    echo ""
    echo "----------------------------------------------"
    echo " Test: ${test_name}"
    echo "----------------------------------------------"

    case "$ENGINE" in
        cromwell)
            if run_test_cromwell "$test_name" "$inputs_json"; then
                PASSED=$((PASSED + 1))
            else
                FAILED=$((FAILED + 1))
            fi
            ;;
        miniwdl)
            if run_test_miniwdl "$test_name" "$inputs_json"; then
                PASSED=$((PASSED + 1))
            else
                FAILED=$((FAILED + 1))
            fi
            ;;
        *)
            echo -e "${RED}Unknown engine: ${ENGINE}. Use 'cromwell' or 'miniwdl'.${NC}"
            exit 1
            ;;
    esac
done

echo ""
echo "=============================================="
echo " Test Results Summary"
echo "=============================================="
echo -e " ${GREEN}Passed:  ${PASSED}${NC}"
echo -e " ${RED}Failed:  ${FAILED}${NC}"
echo "=============================================="

if [[ $FAILED -gt 0 ]]; then
    exit 1
fi
