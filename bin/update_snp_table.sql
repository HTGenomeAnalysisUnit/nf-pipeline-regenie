DROP TABLE IF EXISTS snps;

CREATE TABLE snps (
	"chrom" TEXT NOT NULL,
	"pos" INTEGER NOT NULL,
	"id" TEXT NOT NULL,
	"ref" TEXT,
	"alt" TEXT,
	"rsid" TEXT
);

.mode tabs
.import "snp.table" snps

CREATE INDEX snps_id_idx ON snps(id);
CREATE INDEX snps_rsid_idx ON snps(rsid);