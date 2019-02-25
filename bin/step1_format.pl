open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].out" or die "$!";

while (<In>) {
chomp;
@temp = split (/\s/, $_);
$temp[1] =~/tsplit_out\/(.*)\.2bit/;
$temp[2] =~/qsplit_out\/(.*)\.2bit/;

print Out "$_ &\n";
}

