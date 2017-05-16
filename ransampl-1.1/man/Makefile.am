man_MANS = ransampl.3
noinst_DATA = ransampl.html

%.3 : %.pod
	pod2man -s 3 -c "ransampl manual" -n "$*"      $< > $@

%.html : %.pod
	pod2html --title="ransampl: a C library for random number sampling" --noindex $< > $@
