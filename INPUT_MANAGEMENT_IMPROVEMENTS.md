# Input Management Improvements for nf-pipeline-regenie

## Overview

This document outlines the comprehensive improvements made to the input management system of the nf-pipeline-regenie pipeline. The enhancements implement modern Nextflow best practices using the nf-schema plugin for robust parameter validation, enhanced logging, and institutional profile support.

## Key Improvements

### 1. nf-schema Plugin Integration

**What changed:**
- Replaced manual parameter validation loops with the modern nf-schema plugin
- Updated JSON schema from draft-07 to draft 2020-12 format
- Added comprehensive parameter validation and summary generation

**Benefits:**
- Automatic parameter validation using JSON schema definitions
- Enhanced parameter summaries with better formatting
- Better error messages for invalid parameters
- Standardized validation approach following nf-core best practices

**Configuration added:**
```groovy
plugins {
    id 'nf-schema@2.3.0'
}

validation {
    parametersSchema = 'nextflow_schema.json'
    ignoreParams = []
    defaultIgnoreParams = ['genomes']
}
```

### 2. Enhanced Logging and Help System

**What changed:**
- Implemented comprehensive parameter logging using `paramsSummaryLog()`
- Added detailed pipeline information display
- Added help functionality with `--help` parameter
- Enhanced user experience with better formatted output

**Features:**
- Automatic parameter summary at pipeline start
- Pipeline version, git information, and execution details
- Help text generation from schema definitions
- Professional-looking log output

**Usage:**
```bash
# Show help message
nextflow run nf-pipeline-regenie --help

# The pipeline automatically displays:
# - Parameter summary
# - Pipeline information
# - Git repository details
# - Execution environment info
```

### 3. Institutional Profile Support

**What changed:**
- Added support for nf-core institutional configurations
- Enabled use of community-maintained compute environment profiles
- Added custom config loading from nf-core/configs repository

**Benefits:**
- Easy integration with existing institutional compute environments
- Access to 150+ pre-configured institutional profiles
- Simplified deployment across different HPC systems
- Standardized resource management

**Configuration added:**
```groovy
params {
    // Custom config version and base for institutional profiles
    custom_config_version      = 'master'
    custom_config_base         = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
    
    // Profile information for institutional configurations
    config_profile_description = null
    config_profile_contact     = null
    config_profile_url         = null
}

// Load institutional configs from nf-core/configs repository
includeConfig !System.getenv('NXF_OFFLINE') && params.custom_config_base ? "${params.custom_config_base}/nfcore_custom.config" : "/dev/null"
```

**Usage examples:**
```bash
# Use UPPMAX institutional profile
nextflow run nf-pipeline-regenie -profile uppmax --project myproject

# Use custom institutional config
nextflow run nf-pipeline-regenie \
  --custom_config_base 'https://raw.githubusercontent.com/myinstitution/configs/main' \
  -profile myprofile
```

### 4. Modern Resource Management

**What changed:**
- Replaced deprecated `check_max()` function with modern `resourceLimits`
- Updated all process resource definitions in `conf/base.config`
- Implemented Nextflow's built-in resource limit enforcement

**Benefits:**
- Better resource management and job scheduling
- Automatic resource limit enforcement by Nextflow
- More predictable behavior on different compute environments
- Future-proof resource handling

**Configuration:**
```groovy
process {
    // Resource limits to prevent jobs from requesting more than available
    resourceLimits = [
        cpus: 128,
        memory: 512.GB,
        time: 7.d
    ]
    
    // Process resources now use simple expressions
    cpus   = { 1    * task.attempt }
    memory = { 4.GB * task.attempt }
    time   = { 1.h  * task.attempt }
}
```

## Usage Examples

### Basic Usage with Enhanced Logging
```bash
nextflow run nf-pipeline-regenie \
  --project my_gwas_study \
  --genotypes_build GRCh38 \
  --genotypes_array /path/to/genotypes \
  --phenotypes_filename /path/to/phenotypes.txt
```

### Using Institutional Profiles
```bash
# Using a specific institutional profile
nextflow run nf-pipeline-regenie \
  -profile crick \
  --project my_gwas_study \
  --genotypes_build GRCh38

# Using multiple profiles
nextflow run nf-pipeline-regenie \
  -profile docker,uppmax \
  --project my_gwas_study
```

### Getting Help
```bash
# Display comprehensive help
nextflow run nf-pipeline-regenie --help

# Help includes:
# - Parameter descriptions from schema
# - Required vs optional parameters
# - Parameter formats and patterns
# - Usage examples
```

### Parameter Validation Examples

The pipeline now provides detailed validation messages:

```bash
# Invalid project name (contains spaces)
nextflow run nf-pipeline-regenie --project "my project"
# ERROR: Parameter validation failed!
# - project: Project name must not contain spaces

# Missing required parameter
nextflow run nf-pipeline-regenie --project myproject
# ERROR: The following required parameters are missing:
# - genotypes_build: Genome build must be specified
```

## Benefits Summary

1. **Improved User Experience:**
   - Clear parameter validation with helpful error messages
   - Comprehensive help system
   - Professional logging output
   - Better documentation through schema

2. **Enhanced Reliability:**
   - Robust parameter validation
   - Modern resource management
   - Standardized configuration approach
   - Better error handling

3. **Institutional Integration:**
   - Easy deployment on institutional clusters
   - Access to community-maintained profiles
   - Standardized configuration sharing
   - Simplified HPC integration

4. **Developer Benefits:**
   - Modern Nextflow best practices
   - Maintainable code structure
   - Future-proof architecture
   - Community-standard approaches

## Migration Notes

For existing users of the pipeline:

1. **No breaking changes:** All existing parameter names and functionality remain the same
2. **Enhanced validation:** Some previously accepted invalid parameters may now be caught and rejected
3. **Better logging:** More detailed output provides better visibility into pipeline execution
4. **New features:** Help system and institutional profiles are opt-in enhancements

## Technical Implementation Details

### Files Modified:
- `main.nf`: Added nf-schema plugin integration and enhanced logging
- `nextflow.config`: Added plugin configuration, institutional profile support, and removed deprecated functions
- `conf/base.config`: Updated resource management to use resourceLimits
- `nextflow_schema.json`: Updated to JSON Schema draft 2020-12 format

### Key Functions Added:
- `paramsSummaryLog()`: Comprehensive parameter summary
- `paramsSummaryMap()`: Parameter data for programmatic use
- `paramsHelp()`: Dynamic help text generation

### Dependencies:
- nf-schema plugin v2.3.0 or higher
- Nextflow version compatible with resourceLimits (21.10+)
- Internet access for institutional profile loading (optional)

This implementation represents a significant improvement in pipeline usability, maintainability, and integration capabilities while maintaining full backward compatibility with existing workflows.