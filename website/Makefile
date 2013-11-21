
FILES=docs index license perl_module specs links htslib

install:
	    @for i in $(FILES); do echo "php src/index.php $$i > $$i.html"; php src/index.php $$i > $$i.html; done;
		zcat ../cpp/vcftools.1.gz | groff -mandoc -Thtml > vcftools.1.html
