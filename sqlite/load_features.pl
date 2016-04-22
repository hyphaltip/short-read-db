#!/usr/bin/perl -w

=head1 NAME 

load_features  - load genomic annotation features into SQLlite database

=head1 SYNOPSIS

 load_features --dbdir database_location --db database_name fastafiles gffiles 

=head1 DESCRIPTION

Comand line options:

=head AUTHOR

Jason Stajich, jason_at_bioperl.org

=cut
    
use strict;
use Getopt::Long;
use File::Spec;
use DBI qw(:sql_types);
use Bio::SeqIO;

my %unzip = ('bz2' => 'bzcat',
	     'gz'  => 'zcat');

my $commit_interval = 200_000; # dbh->commit interval, tune based on disk/memory

my $dbname;
my $dbdir  = '/tmp/sReadDB';
my $seqformat = 'fasta';

my $debug = 0;
my $create = 0;
GetOptions(
	   'v|versbose!'      => \$debug,
	   'db|dbname:s'      => \$dbname, # required
	   'd|dir|dbdir:s'    => \$dbdir,
	   'c|create!'        => \$create,
           );

if( ! defined $dbname ) {
    die("must provide a dbname with --db or --dbname\n");
} 
my $dbargs = { AutoCommit => 0,
	       RaiseError => 1,
	       PrintError => 0};

my $dbidx = File::Spec->catfile($dbdir,"$dbname.db");
mkdir($dbdir) unless -d $dbdir;
if( ! -f $dbidx || -z $dbidx ) {
    $create = 1;
}

my $dbh = DBI->connect("dbi:SQLite:$dbidx","","",$dbargs);

if( $create ) {    
    initialize($dbh);
}

my (@gff,@fasta);

for (@ARGV) {
    if (/\.(fa|fna|fasta|dna|seq|fast)(?:$|\.)/i) {
	push @fasta,$_;
    } else {
	push @gff,$_;
    }
}

my $chrom2id = load_chroms($dbh,\@fasta);
load_features($dbh,$chrom2id,\@gff);

sub initialize {
    $dbh->do('DROP TABLE IF EXISTS chromosome');
    $dbh->do(<<EOF 	    
	     CREATE TABLE chromosome 
	     (
	      chromosome_id integer PRIMARY KEY,
	      name varchar(64) NOT NULL,
	      length integer NOT NULL,
	      description varchar(255) default NULL
	      )
EOF
);	     
    $dbh->do('CREATE UNIQUE INDEX ui_chrom_name ON chromosome (name)');
    $dbh->do('CREATE INDEX i_chrom_length ON chromosome (length)');

    $dbh->do('DROP TABLE IF EXISTS feature');
    $dbh->do(<<EOF 	    
	 CREATE TABLE feature 
	 ( 
	   feature_id integer PRIMARY KEY,
	   chromosome_id integer NOT NULL,
	   load_id varchar(128) NOT NULL,
	   fname   varchar(128) default NULL,
	   fparent varchar(128) default NULL,
	   ftype   varchar(64) NOT NULL,
	   fsource varchar(64) NOT NULL,
	   fstart  integer NOT NULL,
	   fstop   integer NOT NULL,
	   fstrand char(1) NOT NULL
	   )
EOF
);	     

    $dbh->do("CREATE INDEX i_feat_chrom ON feature (chromosome_id)");
    $dbh->do("CREATE INDEX i_feat_start ON feature (fstart)");
    $dbh->do("CREATE INDEX i_feat_stop ON feature (fstop)");
    $dbh->do("CREATE INDEX i_feat_name ON feature (fname)");
    $dbh->do("CREATE UNIQUE INDEX i_id ON feature (load_id)");

    $dbh->do('DROP TABLE IF EXISTS read_location');    
    $dbh->commit;
}


sub load_chromosomes {
    my $dbh = shift;
    my $files = shift;
    my %chrom2id;
    my $sth = $dbh->prepare(<<EOF
			    INSERT INTO chromosome (name,length,description)
			    VALUES (?,?,?)
			    
EOF
			    );
	for my $fasta ( @$files ) {
	    my $in;
	    if( $fasta =~ /\.(gz|bz2)$/ ) {
		$in = Bio::SeqIO->new(-format => $seqformat,
				      -file   => "$unzip{$1} $fasta |");
	    } else {
		$in = Bio::SeqIO->new(-format => $seqformat,
				      -file   => $fasta);
	    }

	    while( my $seq = $in->next_seq ) {
		eval {
		    $sth->execute($seq->display_id,
				  $seq->length,
				  $seq->description || '');
		};
		if( $@ ) {
		    warn($@,"\n") if $@ !~ /unique/;
		}
	    }
	    $dbh->commit;
	    $sth->finish;
	}
    
    # chromosome load debug    
    my $sth = $dbh->prepare("SELECT * from chromosome");
    $sth->execute;
    my ($pid,$name,$len,$desc);
    $sth->bind_columns(\$pid,\$name,\$len,\$desc);
    while( $sth->fetch ) {
	$chrom2id{$name} = $pid;
	print join("\t", $pid,$name,$len,$desc), "\n" if $debug;
    }
    $sth->finish;
    \%chrom2id;
}


sub load_features {
    my $dbh = shift;
    my $chrom2id = shift;
    my $files = shift;
    my $sth = $dbh->prepare(<<EOF
INSERT INTO feature (chromosome_id,load_id,
		     fname,fparent,
		     ftype,fsource,
		     fstart,fstop,fstrand)
	     VALUES (?,?,?,?,?,?,?,?,?)
EOF			    
			    );
    for my $gff ( @$files ) {
	my $fh;
	if( $gff =~ /\.(gz|bz2)$/ ) {
	    open( $fh => "$unzip{$1} $gff |") || die $!;
	} else {
	    open($fh => "< $gff") || die $!;
	}
	while(<$gff>) {
	    next if /^\#/ || /^\s+$/;
	    chomp;
	    my ($seqid,$src,$type,$start,$end,
		$score,
		$strand,
		undef,
		$group) = split(/\t/,$_);
	    my %group = map { split(/=/,$_) } split(/;/,$group);
	    if( $group{'ID'} ) {
	    }
	}
    }
    
}
END {
    $dbh->commit;
    $dbh->disconnect;
}
