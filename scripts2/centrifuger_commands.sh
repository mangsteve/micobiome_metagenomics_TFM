

centrifuger/centrifuger/centrifuger-download -o taxonomy taxonomy

centrifuger/centrifuger/centrifuger-download -o library -d "fungi" refseq > seqid2taxid.map

ls library/*/*.fna.gz > file.list 

centrifuger/centrifuger/centrifuger-build -t 4 --conversion-table seqid2taxid.map \
	--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
	-l file.list -o refseq_abv --build-mem 14G

    