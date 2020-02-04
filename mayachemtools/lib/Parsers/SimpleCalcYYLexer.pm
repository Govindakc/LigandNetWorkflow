package Parsers::SimpleCalcYYLexer;
#
# File: SimpleCalcYYLexer.pm
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2019 Manish Sud. All rights reserved.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#

use strict;
use Carp;
use Exporter;
use Scalar::Util ();
use Parsers::YYLexer;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Parsers::YYLexer Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName, $YYTabFile, @YYLexerTokensSpec);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifySimpleCalcYYLexer';

# Class constructor...
sub new {
  my($Class, $Input) = @_;
  my(@TokensSpec);

  # Initialize object...
  my $This = $Class->SUPER::new($Input,  @YYLexerTokensSpec);
  bless $This, ref($Class) || $Class;
  $This->_InitializeYYLexer();

  return $This;
}

# Initialize object data...
#
sub _InitializeYYLexer {
  my($This) = @_;

  # Setup default YYTabFile containing mapping of token names to numbers...
  $This->SetupYYTabFile($YYTabFile);

  return $This;
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;

  # Setup default token table file...
  $YYTabFile = "Parsers/SimpleCalcParser.tab.ph";

  # Setup default lexer tokens specs...
  @YYLexerTokensSpec = (
			[ 'LETTER', qr/[a-zA-Z]/ ],
			[ 'NUMBER', qr/\d+/ ],
			[ 'SPACE', qr/[ ]*/, sub { my($This, $TokenLabel, $MatchedText) = @_; return ''; } ],
			[ 'NEWLINE', qr/(?:\r\n|\r|\n)/, sub { my($This, $TokenLabel, $MatchedText) = @_;  return "\n"; } ],
			[ 'CHAR', qr/./ ]
		       );
}

# Is it a lexer object?
sub _IsSimpleCalcYYLexer {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Return a string containing information about lexer...
sub StringifySimpleCalcYYLexer {
  my($This) = @_;
  my($SimleCalcYYLexerString);

  $SimleCalcYYLexerString = "SimpleCalcYYLexer: PackageName: $ClassName; " . $This->_GetYYLexerInfoString();

  return $SimleCalcYYLexerString;
}

1;

__END__

=head1 NAME

Parsers::SimpleCalcYYLexer

=head1 SYNOPSIS

use Parsers::SimpleCalcYYLexer;

use Parsers::SimpleCalcYYLexer qw(:all);

=head1 DESCRIPTION

B<SimpleCalcYYLexer> class provides the following methods:

new, YYLex, GetYYLex, StringifySimpleCalcYYLexer

B<Parser::SimpleCalcYYLexer> class is derived from B<Parser::YYLexer> class, which in
turn is derived from base class B<Parser::Lexer> that provides all the underlying
lexer functionality. B<SimpleCalcYYLexer> class is designed to be used with
B<yyparse> code generated by running B<byacc> on a parser defined using
parser definition B<SimpleCalcParser.yy> file.

The parser package and token table files, B<SimpleCalcParser.pm> and B<SimpleCalcParser.tab.ph>,
are automatically generated from parser grammar definition file, B<SimpleCalcParser.yy>, using
byacc available through perl-byacc1.8 modified with perl5-byacc-patches-0.5 for generation
of object oriented parser:

    byacc -l -P -d -b SimpleCalcParser SimpleCalcParser.yy
    mv SimpleCalcParser.tab.pl SimpleCalcParser.pm

B<SimpleCalcYYLexer.pm> class implements a lexer for a simple calculator and is provided
to highlight usasge of B<YYLex> through B<yyparse>.

The default specification of lexer tokens for B<SimpleCalcYYLexer.pm> includes:

    @YYLexerTokensSpec = (
        [ 'LETTER', qr/[a-zA-Z]/ ],
        [ 'NUMBER', qr/\d+/ ],
        [ 'SPACE', qr/[ ]*/,
            sub { my($This, $TokenLabel, $MatchedText) = @_; return ''; }
        ],
        [ 'NEWLINE', qr/(?:\r\n|\r|\n)/,
            sub { my($This, $TokenLabel, $MatchedText) = @_;  return "\n"; }
        ],
        [ 'CHAR', qr/./ ]
    );

The default B<SimpleCalcParser.tab.ph> file containing token identifiers for
B<SimpleCalcParser.yy> includes:

    $NUMBER=257;
    $LETTER=258;

=head2 METHODS

=over 4

=item B<new>

    $SimpleCalcYYLexer = new Parsers::SimpleCalcYYLexer($Input);

Using specified I<Input>, B<new> method generates a new B<SimpleCalcYYLexer>
and returns a reference to newly created B<SimpleCalcYYLexer> object.

Examples:

    # Input string...
    $InputText = "3 + 4 +6\nx=3\ny=5\nx+y\nx+z\n";

    $YYLexer = new Parsers::SimpleCalcYYLexer($InputText);
    $YYLex = $YYLexer->GetYYLex();

    $Debug = 0;
    $CalcParser = new Parsers::SimpleCalcParser($YYLex,
                            \&Parsers::SimpleCalcParser::yyerror, $Debug);
    $Value = $SimpleCalcParser->yyparse();
    print "Value: $Value\n";

    # Input file...
    $InputFile = "Input.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $YYLexer = new Parsers::SimpleCalcYYLexer($InputFile);
    $YYLex = $YYLexer->GetYYLex();

    $CalcParser = new Parsers::SimpleCalcParser($YYLex,
                            \&Parsers::SimpleCalcParser::yyerror);
    $Value = $SimpleCalcParser->yyparse();
    print "Value: $Value\n";

    # Input file iterator...
    $InputFile = "TestSimpleCalcParser.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $InputIterator = sub { return <INPUTFILE>; };
    $YYLexer = new Parsers::SimpleCalcYYLexer($InputIterator);
    $YYLex = $YYLexer->GetYYLex();

    $CalcParser = new Parsers::SimpleCalcParser($YYLex,
                            \&Parsers::SimpleCalcParser::yyerror);
    $Value = $SimpleCalcParser->yyparse();
    print "Value: $Value\n";

=item B<StringifySimpleCalcYYLexer>

    $YYLexerString = $YYLexer->StringifySimpleCalcYYLexer();

Returns a string containing information about I<YYLexer> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Lexer.pm, YYLexer.pm, SimpleCalcParser.yy

=head1 COPYRIGHT

Copyright (C) 2019 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut