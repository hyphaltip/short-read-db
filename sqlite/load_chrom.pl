#!/usr/bin/perl -w
use strict;

=head1 NAME 

load_chrom  - load chromosome info into SQLlite database

=head1 SYNOPSIS

 load_features --dbdir database_location --db database_name --fasta chromosomes

=head1 DESCRIPTION

Comand line options:

=head AUTHOR

Jason Stajich, jason_at_bioperl.org

=cut
    
use strict;
use Getopt::Long;
use File::Spec;
use DBI qw(:sql_types);

my %unzip = ('bz2' => 'bzcat',
	     'gz'  => 'zcat');

my $commit_interval = 200_000; # dbh->commit interval, tune based on disk/memory

my $dbname;
my $dbdir  = '/tmp/sReadDB';

my $debug  = 0;
my $create = 0;
GetOptions(
	   'v|versbose!'      => \$debug,
	   'db|dbname:s'      => \$dbname, # required
	   'd|dir|dbdir:s'    => \$dbdir,
	   'i|input|gff:s'    => \@input,
	   'c|create!'        => \$create,
           );

if( ! defined $dbname ) {
    die("must provide a dbname with --db or --dbname\n");
} 
my $dbargs = { AutoCommit => 0,
	       PrintError => 1};

if( ! @input ) {
    warn("no input files\n");
    exit;
}

mkdir($dbdir) unless -d $dbdir;

my $dbidx = File::Spec->catfile($dbdir,"$dbname.db");

$dbh = DBI->connect("dbi:SQLite:$dbidx","","",$dbargs);

if( $create ) {    
    $dbh->do('DROP TABLE chromosome') if $clean;
    $dbh->do('DROP TABLE feature') if $clean;
    $dbh->do('DROP TABLE read_location') if $clean;
    $dbh->do(<<EOF
	 CREATE TABLE chromosome 
	 ( 
	   feature_id integer PRIMARY KEY,
	   chromosome_id integer NOT NULL,
	   load_id varchar(128) NOT NULL,
	   fname   varchar(128) default NULL,
	   fparent varchar(128) default NULL,
	   ftype   varchar(64) NOT NULL,
	   fstart  integer NOT NULL,
	   fstop   integer NOT NULL,
	   fstrand integer NOT NULL
	   ));
$dbh->do("CREATE INDEX i_chrom ON feature (chromosome_id)");
$dbh->do("CREATE INDEX i_start ON feature (fstart)");
$dbh->do("CREATE INDEX i_stop ON feature (fstop)");
$dbh->do("CREATE INDEX i_name ON feature (fname)");
$dbh->do("CREATE UNIQUE INDEX i_id ON feature (load_id)");

}

