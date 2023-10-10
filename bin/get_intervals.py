import sqlite3

# Connect to the database
conn = sqlite3.connect('your_database.db')
c = conn.cursor()

# Define the interval size
interval_size = 100

# Get the distinct chromosomes in the Variant table
c.execute('SELECT DISTINCT chromosome FROM Variant')
chromosomes = [row[0] for row in c.fetchall()]

# Loop over the chromosomes and generate intervals
intervals = []
for chromosome in chromosomes:
	c.execute(f'SELECT position FROM Variant WHERE chromosome = "{chromosome}" ORDER BY position ASC')
	positions = [row[0] for row in c.fetchall()]
	for i in range(0, len(positions), interval_size):
		start = positions[i]
		end = positions[min(i + interval_size - 1, len(positions) - 1)]
		interval_name = f'{chromosome}:{start}-{end}'
		intervals.append(interval_name)

# Save the intervals to a file
with open('intervals.txt', 'w') as f:
	f.write('\n'.join(intervals))

# Close the database connection
conn.close()