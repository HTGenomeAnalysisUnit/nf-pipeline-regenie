.mode csv
.output intervals.txt

WITH
	-- Get the distinct chromosomes in the Variant table
	chromosomes AS (
		SELECT DISTINCT chromosome FROM Variant
	),
	-- Get the positions of the variants for each chromosome
	positions AS (
		SELECT chromosome, position FROM Variant
	),
	-- Generate intervals of size %CHUNK_SIZE% for each chromosome
	intervals AS (
		SELECT
			chromosome,
			MIN(position) AS start,
			MAX(position) AS end
		FROM (
			SELECT
				chromosome,
				position,
				(ROW_NUMBER() OVER (PARTITION BY chromosome ORDER BY position) - 1) / 100 AS interval
			FROM positions
		)
		GROUP BY chromosome, interval
	)
-- Generate the interval names and output to file
SELECT
	chromosome || ':' || start || '-' || end AS interval_name
FROM intervals
ORDER BY chromosome, start;

.quit