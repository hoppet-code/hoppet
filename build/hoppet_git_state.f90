module hoppet_git_state
  implicit none  

  character(len=*), parameter :: HOPPET_GIT_STATUS_UNO = " M .github/workflows/build-wheels.yml;"//&
     " M pyinterface/CMakeLists.txt;"//&
     " M pyinterface/setup.py;"//&
     " M pyproject.toml"
  character(len=*), parameter :: HOPPET_GIT_AUTHOR_NAME = "copilot-swe-agent[bot]"
  character(len=*), parameter :: HOPPET_GIT_AUTHOR_EMAIL = "198982749+Copilot@users.noreply.github.com"
  character(len=*), parameter :: HOPPET_GIT_HEAD_SHA1 = "6ab1d16b35047769df6d46c8475c75de5516c134"
  character(len=*), parameter :: HOPPET_GIT_COMMIT_DATE_ISO8601 = "2025-09-25 14:03:41 +0000"
  !character(len=*), parameter :: HOPPET_GIT_COMMIT_SUBJECT = "Initial plan"
  character(len=*), parameter :: HOPPET_GIT_DESCRIBE = "6ab1d16"
  !logical, parameter :: GIT_RETRIEVED_STATE = true
  !logical, parameter :: GIT_IS_DIRTY = true

contains

  subroutine hoppet_print_git_state(idev,prefix)
    integer, intent(in) :: idev
    character(len=*), intent(in), optional :: prefix

    if (present(prefix)) then
        write(idev,'(a)',advance='no') trim(prefix)//" "
    end if
    write(idev,'(a)') "hoppet compiled git info: "//trim(HOPPET_GIT_HEAD_SHA1)&
       &              // " , " // HOPPET_GIT_COMMIT_DATE_ISO8601 // ", status: " &
       &              // trim(HOPPET_GIT_STATUS_UNO)
  end subroutine hoppet_print_git_state

end module hoppet_git_state

! bool GitMetadata::Populated() {
!     return true;
! }
! bool GitMetadata::AnyUncommittedChanges() {
!     return true;
! }
! std::string GitMetadata::StatusUNO() {
!     return " M .github/workflows/build-wheels.yml; M pyinterface/CMakeLists.txt; M pyinterface/setup.py; M pyproject.toml";
! }
! std::string GitMetadata::AuthorName() {
!     return "copilot-swe-agent[bot]";
! }
! std::string GitMetadata::AuthorEmail() {
!     return "198982749+Copilot@users.noreply.github.com";
! }
! std::string GitMetadata::CommitSHA1() {
!     return "6ab1d16b35047769df6d46c8475c75de5516c134";
! }
! std::string GitMetadata::CommitDate() {
!     return "2025-09-25 14:03:41 +0000";
! }
! std::string GitMetadata::CommitSubject() {
!     return "Initial plan";
! }
! std::string GitMetadata::Describe() {
!     return "6ab1d16";
! }
