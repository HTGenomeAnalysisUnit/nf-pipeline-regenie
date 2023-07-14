import gwaslab as gl
gl.options.set_option('data_directory', '/gwaslab_data/')

gl.download_ref('1kg_eur_hg19')
gl.download_ref('1kg_eur_hg38')
gl.download_ref('1kg_eur_hg19_tbi')
gl.download_ref('recombination_hg19')
gl.download_ref('recombination_hg38')
gl.download_ref('ensembl_hg19_gtf')
gl.download_ref('ensembl_hg38_gtf')
