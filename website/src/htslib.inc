
<div class="item">
<h1>The bcftools/htslib VCF commands</h1>
<p> <a href="https://github.com/samtools/htslib">HTSlib</a> is a C library for
high-throughput sequencing data formats. It is designed for speed and works with both VCF and 
<A href="http://samtools.github.io/hts-specs/">BCFv2</A>.
</p>

<h3>Download and installation</h3>
<p>
The library is hosted on github. It can be downloaded and compiled the usual way.
The <span class="cmd">clone</span> command is run only once, the <span class="cmd">pull</span>
command is run whenever the latest snapshots from github is needed.
Please see the <a href="https://github.com/samtools/bcftools">bcftools github
page</a> for the up-to-date version of the clone command. The software is under heavy
development and the option <span class="cmd">--branch</span> may be required.
</p>
<p class="codebox">
     git clone [<i>--branch=name</i>] git://github.com/samtools/htslib.git htslib
<br> git clone git://github.com/samtools/bcftools.git bcftools 
<br> cd htslib; git pull; cd ..
<br> cd bcftools; git pull; cd ..
<br>
<br> # Compile
<br> cd bcftools; make; make test
<br>
<br> # Run
<br> ./bcftools stats file.vcf.gz
<br> 
</p>

<h3>The tools</h3>
<ul class="faq">
<li> <a href="#annotate">bcftools annotate</a>
<li> <a href="#call">bcftools call</a>
<li> <a href="#filter">bcftools filter</a>
<li> <a href="#gtcheck">bcftools gtcheck</a>
<li> <a href="#isec">bcftools isec</a>
<li> <a href="#merge">bcftools merge</a>
<li> <a href="#norm">bcftools norm</a>
<li> <a href="#query">bcftools query</a>
<li> <a href="#stats">bcftools stats</a>
<li> <a href="#view">bcftools view</a>
</ul>

</div>

<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
<script type="text/javascript">
jQuery(document).ready(function() {
        jQuery(".usageText").hide();
        jQuery(".usageToggle").click(function() { jQuery(this).next(".usageText").slideToggle(100); });
    });
</script>

<div class="item">
<a name="annotate" class="Q">bcftools annotate</a>
<p> Adds or removes annotations, support for user-written plugins.
</p>
<p> Fast alternative to <A href="perl_module.html#vcf-annotate">vcf-annotate</A>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools annotate 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="call" class="Q">bcftools call</a>
<p> Formerly known as <span class="cmd">bcftools view</span>, this is the successor of the popular caller from the <span class="cmd">samtools</span> package with extended capabilities. 
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools call 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="filter" class="Q">bcftools filter</a>
<p> Powerful fixed-threshold filtering, accepts boolean and arithmetic expressions.
<br>
See also the <A href="view">bcftools view</A> below.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools filter 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="gtcheck" class="Q">bcftools gtcheck</a>
<p> A tool for detecting sample swaps and contamination
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools gtcheck 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="isec" class="Q">bcftools isec</a>
<p> Fast alternative to <A href="perl_module.html#vcf-isec">vcf-isec</A>
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools isec 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="merge" class="Q">bcftools merge</a>
<p> Fast alternative to <A href="perl_module.html#vcf-merge">vcf-merge</A> with extended capabilities and correct handling of Number=A,G,R INFO fields.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools merge 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="norm" class="Q">bcftools norm</a>
<p> Left-align and normalize indels to the shortest possible representation.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools norm 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="query" class="Q">bcftools query</a>
<p> Fast alternative to <A href="perl_module.html#vcf-query">vcf-query</A>
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools query 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="stats" class="Q">bcftools stats</a> (former vcfcheck)
<p> Extract stats from a VCF/BCF file or compare two VCF/BCF files. The resulting text file can be plotted using 
<span class="cmd">plot-vcfstats</span>.
</p>
<p class="codebox">
bcftools stats file.vcf.gz &gt; file.vchk 
<br> plot-vcfstats file.vchk -p plots/</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools stats 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>

<div class="item">
<a name="view" class="Q">bcftools view</a>
<p> This versatile tool can be used for subsetting by sample, position and even flexible fixed-threshold filtering.
</p>
<div class="usageBox"><span class="usageToggle">(Read more)</span><pre class="usageText"><?php unset($out);exec("bcftools view 2>&1",$out); echo htmlspecialchars(implode("\n",$out)); ?></pre></div>
</div>


