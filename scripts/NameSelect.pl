#!/usr/bin/perl -w

# perl script which parses STDIN (for now) and  produces
# as output a file which allows conversion between codes and
# strings
#
# Author: Gavin P. Salam, 2002.

# GPS mods 31/07/2003 & 01/08:
#
#   . added option "-prefix PrefixString" which causes all module
#     and public subprogram names to be prefixed with PrefixString
#
#   . comments at beginning of module and subprograms now use !! so
#     as to be picked up by f90doc
#
#   . added usage string
#
#   . added optional arg longname to NameOfCode. The longname is taken
#     (a) as any string following !! on the same line as the constant or 
#     if such a string is not present (b) any string following !! on the 
#     previous line
#

# modules that will be used
%modules = ();
%longnames = ();
@consts = ();

# remember that these can be modified by prefix list...
$ThisModuleName="NameSelect";
$CodeOfName="CodeOfName";
$NameOfCode="NameOfCode";
$code_val_opt="code_val_opt";


if ($#ARGV < 0) {
  print "Usage:\n";
  print "  NameSelect.pl [-prefix PREFIX] file1.f90 [file2.f90 [...]]\n";
  exit 0;
}

while ($arg = shift @ARGV) {
  if ($arg =~ /^-/) {
    if ($arg eq '-prefix') {
      $prefix = shift @ARGV;
      $ThisModuleName = $prefix.$ThisModuleName;
      $CodeOfName     =	$prefix.$CodeOfName    ;
      $NameOfCode     =	$prefix.$NameOfCode    ;
      $code_val_opt   =	$prefix.$code_val_opt  ;
      next;
    } else {die "Unrecognized command-line argument $arg"}
  }
  $InputFile = $arg;
  open INPUT, "<$InputFile" || die "Could not open $InputFile\n";
  $last_line = '';
  while ($line = <INPUT>) {
    chomp($line);
    #-- record start of a new module
    if ($line =~ /^[[:space:]]*module/i) {
      ($modname = $line) =~ 
	s/^[[:space:]]*module[[:space:]]*([^![:space:]]*).*/$1/i;
      #print STDERR "Starting module: $modname\n";
    }
    
    if ($line =~ /^[[:space:]]*integer[[:space:],]*public[[:space:],]*parameter/i || $line =~ /^[[:space:]]*integer[[:space:],]*parameter[[:space:],]*public/i)
      {
	$modules{$modname} = 1;
	($const_name = $line) =~ s/[^!]*::[[:space:]]*([^=![:space:]]*).*/$1/i;
	#print STDERR "constant is: $const_name\n";
	push @consts, $const_name;
	# allow possibility of adding a long description of the const
	$longname = '';
	if ($line =~ /!!.+/) {
	  ($longname = $line) =~ s/[^!]*!![[:space:]]*//;
	} elsif ($last_line =~ /^[[:space:]]*!!.+/) {
	  ($longname = $last_line) =~ s/^[[:space:]]*!![[:space:]]*//;
	}
	$longname =~ s/[[:space:]]*$//;
	$longname =~ s/\'//g;
	$longnames{$const_name} = $longname if $longname;
      }
    $last_line = $line;
  }
  close(INPUT);
}

print STDERR "Automatic generation of $ThisModuleName\n";

#print $#consts,"\n";
if ($#consts < 0) {
  print STDERR "No constants were found among the input files\n";
  exit}


open (MODULE, "> $ThisModuleName.f90") || die "Could not open output module";

print MODULE "
!!
!! Module generated automatically by NameSelect.pl
!! Provides conversion from string codes to integer codes
!!
module $ThisModuleName\n";
foreach $module_name (keys %modules) {
  print MODULE "  use $module_name\n";
}

print MODULE "  
  implicit none
  private

  public :: $CodeOfName
  public :: $NameOfCode
  public :: $code_val_opt

  integer, parameter :: maxlen_name = 40
  integer, parameter :: maxlen_longname = 120
  integer, parameter, public :: ${ThisModuleName}_maxlen_name = maxlen_name
  integer, parameter, public :: ${ThisModuleName}_maxlen_longname = maxlen_longname

contains

  !!
  !! Returns the code associated with 'name'; if status dummy arg
  !! is present a non-zero result indicates failure to find the name.
  !! Otherwise failure results in a hard error.
  !!
  function $CodeOfName(name, status) result(code)
    character(len=*), intent(in)   :: name
    integer                        :: code
    integer, optional, intent(out) :: status

    if (present(status)) status = 0

    ! now follows code generated specifically for above used modules
    select case(trim(name))
";

foreach $const_name (@consts) {
  print MODULE 
"    case('$const_name')
      code = $const_name
";}

print MODULE
"    case default
      if (present(status)) then
        status = -1
      else
        write(0,*) 'ERROR. $CodeOfName: unrecognized name \"'//name//'\"'
        stop
      endif
    end select
  end function $CodeOfName


  ! standard code (independent of used modules)
  !!
  !! looks for the command-line argument 'option' and if it is
  !! present its value (optionally prefixed with 'prefix') is
  !! fed to $CodeOfName
  !!
  function $code_val_opt(option, default, prefix) result(code)
    use sub_defs_io
    character(len=*),           intent(in) :: option
    integer, optional,          intent(in) :: default
    character(len=*), optional, intent(in) :: prefix
    integer                                :: code
    !---------------------------------------
    character(len=maxlen_name) :: opt_val

    if (log_val_opt(option)) then
       opt_val = string_val_opt(option)
       if (present(prefix)) then
         code = $CodeOfName(prefix//trim(opt_val))
       else
         code = $CodeOfName(trim(opt_val))
       end if
    else
       if (present(default)) then
          code = default
       else
          write(0,*) 'Error in $code_val_opt: command-line option '&
               &//option//' absent and no default provided'
       end if
    end if
  end function $code_val_opt

  !!
  !! Returns the name of the given integer code which has
  !! the (optional) specified prefix
  !!
  function $NameOfCode(code,prefix,longname) result(name)
    integer,          intent(in)            :: code
    character(len=*), intent(in),  optional :: prefix
    character(len=*), intent(out), optional :: longname
    character(len=maxlen_name)              :: name
    !----------------------------------------------
    integer :: nocc

    nocc = 0
    if (present(longname)) longname = ''

    ! code specific to used modules starts here";

foreach $const_name (@consts) {
  print MODULE 
"
    if (PrefixMatches('$const_name',prefix) .and. code == $const_name) then
       nocc = nocc + 1; name = '$const_name'
";
  if (exists $longnames{$const_name}) {
    $longname = $longnames{$const_name};
    #print STDERR "Writing $longname\n";
    print MODULE
"       if (present(longname)) longname = '$longname'
";
  }
  print MODULE
"    end if
";}
 
print MODULE "
    ! common-code resumes
    if (nocc == 0) then
       if (present(prefix)) then
          write(0,*) 'Error in $NameOfCode: could not find code',&
               &code,' with prefix \"'//prefix//'\"'
       else
          write(0,*) 'Error in $NameOfCode: could not find code',&
               &code,' (without prefix)'
       end if
       stop
    end if

    if (nocc > 1) then
       if (present(prefix)) then
          write(0,*) 'Error in $NameOfCode: several meanings for code',&
               &code,' with prefix \"'//prefix//'\"'
       else
          write(0,*) 'Error in $NameOfCode: several meanings for code',&
               &code,' (without prefix)'
       end if
       stop
    end if

  end function $NameOfCode



  !!
  !! for establishing whether a prefix matches a given string
  !!
  logical function PrefixMatches(string,prefix)
    character(len=*), intent(in)           :: string
    character(len=*), intent(in), optional :: prefix

    if (present(prefix)) then
       PrefixMatches = (index(string,prefix) == 1)
    else
       PrefixMatches = .true.
    end if
  end function PrefixMatches


end module $ThisModuleName
";

close(MODULE);

