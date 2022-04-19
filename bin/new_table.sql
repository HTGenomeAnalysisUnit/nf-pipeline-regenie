DROP TABLE IF EXISTS @table_name@;

CREATE TABLE @table_name@ (
    "RUN_ID" TEXT,
    "PHENO" TEXT,
    "CHROM" TEXT,
    "GENPOS" INTEGER,
    "ID" TEXT,
    "ALLELE0" TEXT,
    "ALLELE1" TEXT,
    "A1FREQ" REAL,
    "INFO" REAL,
    "N" INTEGER,
    "TEST" TEXT,
    "EFFECT" REAL,
    "SE"    REAL,
    "CHISQ" REAL,
    "LOG10P"    REAL,
    "EXTRA" TEXT
);

.mode tabs
.import "| tail -n+2 @filename@" @table_name@

CREATE INDEX @table_name@_snpid ON @table_name@ (ID);
CREATE INDEX @table_name@_pos ON @table_name@ (CHROM, GENPOS);
CREATE INDEX @table_name@_pval ON @table_name@ (LOG10P);