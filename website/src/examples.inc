<body>
<div id="wrap">
  <div id="top">
    <h2> <a href="index.html"><strong>VCFtools</strong></a></h2>
    <div id="menu">
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><a href="http://sourceforge.net/projects/vcftools/">Sourceforge page</a></li>	
        <li><a href="examples.html" class="current">Examples &amp; Documentation</a></li>
        <li><a href="downloads.html">Downloads</a></li>
      </ul>
    </div>
  </div>
  <div id="content">
    <div id="left">
    	<h2>Getting Started with VCFtools</h2>
		<p>First, download the program by one of the methods found on the <a href="downloads.html">downloads</a> page.
		Open a terminal and change the working directory to the location of the file download.</p>
		<p>If the program was downloaded through sourceforge, it must be uncompressed using the following command (replacing it with current version number):</p>
		<p class="codebox">tar -xvf vcfools.0.X.XX.tar.gz</p>
		<p>Now change to the vcftools directory and compile the program. If the program was instead downloaded through subversion, change to the specified target directory.</p>
		<p class="codebox">cd vcfools/
		<br>make</p>
		<p>This will compile both the PERL API and the C++ executable. By default, the compiled programs will now be in the <b>vcftools/bin/</b> directory.
		Some common compiler errors may appear, such as not being able to locate <b>Vcf.pm</b> or not being able to locate <b>zlib</b>.
		Please make sure zlib is installed and that the <b>PERL5LIB</b> path is set correctly.</p>
    	<h2>Specific examples by module</h2>
      <h2><a href="perl_module.html">Usage Examples for the PERL API</a></h2>
      <h2><a href="documentation.html">Usage Examples for the Binary Executable</a></h2>
    </div>
     <div id="right">
      <div class="box">
      	<h2>Binary executable documentation</h2>
      		<p>The manual page for <a href="man.html">v0.1.12</a></p>
      	<h2>HTSlib</h2>
      		<p>Documentation for the separate HTSlib package is <a href="htslib.html">here</a>.</p>
      </div>
    </div>
    <div id="clear"></div>
  </div>
  <div id="footer">
    <p>Copyright 2014 &copy; VCFtools</p>
  </div>
</div>
</body>
</html>
