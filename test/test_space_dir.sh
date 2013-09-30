#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

config_dir=$(readlink -f "${BASEDIR}/../config")
space_dir="${BASEDIR}/space dir"
space_dir_touched_file="${space_dir}/touched"
space_dir_makefile="${space_dir}/Makefile"
space_dir_makefile_tmpl="${space_dir}/Makefile.in"
error=0

cleanup() {
    if [ -e "${space_dir_makefile}" ]; then
        cd "${space_dir}" && make clean && cd -
    fi
    rm -rf "${space_dir_touched_file}"
    rm -rf "${space_dir_makefile}"
}

create_makefile() {
    sed -e "s|@CONFIG_DIR@|${config_dir}|" \
        < "${space_dir_makefile_tmpl}" > "${space_dir_makefile}"
}

cleanup
create_makefile

cd "${space_dir}" && make && cd -

if [ -e "${space_dir_touched_file}" ]; then
    echo "Info: file found: ${space_dir_touched_file}"
else
    echo "Error: file not found: ${space_dir_touched_file}"
    error=1 
fi

if [ -e "${space_dir_touched_file}_1" ]; then
    echo "Info: file found: ${space_dir_touched_file}_1"
else
    echo "Warning: file not found: ${space_dir_touched_file}_1"
fi
if [ -e "${space_dir_touched_file}_2" ]; then
    echo "Info: file found: ${space_dir_touched_file}_2"
else
    echo "Warning: file not found: ${space_dir_touched_file}_2"
fi

cleanup

exit ${error}
