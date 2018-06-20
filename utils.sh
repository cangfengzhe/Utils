# set  color
RED="$(tput setaf 1)"
GREEN="$(tput setaf 2)"
YELLOW="$(tput setaf 3)"
BLUE="$(tput setaf 4)"
BOLD="$(tput bold)"
NORMAL="$(tput sgr0)"

check() {
        eval "${1}"
        RESULT=$?
        if [ "${RESULT}" -ne 0 ]
        then
                message "{\"subject\":\"${ME}\",\"input\":${JSON},\"output\":{\"error\":{\"command\":\"${1}\",\"code\":${RESULT}}}}"
                elapsed
                exit ${RESULT}
        fi
}

checkfile() {
        if [ ! -e "${1}" ]
        then
                message "{\"subject\":\"${ME}\",\"input\":${JSON},\"output\":{\"error\":\"Cannot find file ${1}\"}}"
                elapsed
                exit 1
        fi
}
