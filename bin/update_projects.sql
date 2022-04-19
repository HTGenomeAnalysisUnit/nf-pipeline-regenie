CREATE TABLE IF NOT EXISTS projects (
	"project_id" TEXT NOT NULL,
	"phenotype_id" TEXT NOT NULL,
	"phenotype_type" TEXT NOT NULL,
	"covariates" TEXT,
	"n_samples" INTEGER,
	"pipeline_version" TEXT,
	"created_user" TEXT,
	"created_data" DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS projects_pheno ON projects(phenotype_id);
CREATE INDEX IF NOT EXISTS projects_id ON projects(project_id);

INSERT INTO projects VALUES