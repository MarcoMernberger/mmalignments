from bioservices import BioMart

# Initialize BioMart connection
bm = BioMart(host="www.ensembl.org")

# Configure query for Mouse genes
bm.new_query()
bm.add_dataset_to_xml("mmusculus_gene_ensembl")

# Add the Ensembl IDs you want to convert
my_ids = ["ENSMUSG00000059552"]
bm.add_filter_to_xml("ensembl_gene_id", ",".join(my_ids))

# Request corresponding RefSeq mapping attributes
bm.add_attribute_to_xml("ensembl_gene_id")
bm.add_attribute_to_xml("ensembl_transcript_id")
bm.add_attribute_to_xml("ensembl_peptide_id")
bm.add_attribute_to_xml("refseq_mrna")
bm.add_attribute_to_xml("refseq_peptide")

# Execute and parse results
xml_query = bm.get_xml()
result_df = bm.query(xml_query)
print(result_df)
