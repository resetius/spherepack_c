my $text = "";

while ($name = glob("*.c")) {
	open FILE, '<', $name;
	while (<FILE>) {
		$text .= $_;
	}
	close FILE;
}

my @funcs = $text =~ m/^(\/\* Subroutine \*\/.*?\))/gms;

print "
#ifndef SPHEREPACK_H
#define SPHEREPACK_H

#include <f2c.h>

#ifdef __cplusplus
extern \"C\" {
#endif
";

foreach my $i (@funcs) {
	print "$i;\n\n";
}

print "

#ifdef __cplusplus
}
#endif

#endif /*SPHEREPACK_H*/
";

