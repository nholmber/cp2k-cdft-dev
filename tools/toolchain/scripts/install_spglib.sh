#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"
spglib_ver=${spglib_ver:-1.10.4}

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_spglib=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_spglib" ] && rm "${BUILDDIR}/setup_spglib"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_spglib" in
    __INSTALL__)
        echo "==================== Installing spglib ===================="
        pkg_install_dir="${INSTALLDIR}/spglib-${spglib_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "spglib-${spglib_ver} is already installed, skipping it."
        else
            if [ -f spglib-${spglib_ver}.tar.gz ] ; then
                echo "spglib-${spglib_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/spglib-${spglib_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d spglib-${spglib_ver} ] && rm -rf spglib-${spglib_ver}
            tar -xzf spglib-${spglib_ver}.tar.gz
            cd spglib-${spglib_ver}
            mkdir build
            cd build
            cmake -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" -DCMAKE_BUILD_TYPE=Release .. > configure.log 2>&1
            make > make.log 2>&1
            make install >> make.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi

        SPGLIB_CFLAGS="-I${pkg_install_dir}/include"
        SPGLIB_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding spglib from system paths ===================="
        check_command pkg-config --modversion spglib
        add_include_from_paths SPGLIB_CFLAGS "spglib.h" $INCLUDE_PATHS
        add_lib_from_paths SPGLIB_LDFLAGS "libspglib.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking spglib to user paths ===================="
        pkg_install_dir="$with_spglib"
        check_dir "$pkg_install_dir/lib"
        check_dir "$pkg_install_dir/lib64"
        check_dir "$pkg_install_dir/include"
        SPGLIB_CFLAGS="-I'${pkg_install_dir}/include'"
        SPGLIB_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_spglib" != "__DONTUSE__" ] ; then
    SPGLIB_LIBS="-lsymspg"
    if [ "$with_spglib" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_spglib"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    fi
       cat <<EOF > "${BUILDDIR}/setup_spglib"
export SPGLIB_CFLAGS="-I$pkg_install_dir/include ${SPGLIB_CFLAGS}"
export SPGLIB_LDFLAGS="${SPGLIB_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__SPGLIB"
export CP_CFLAGS="\${CP_CFLAGS} ${SPGLIB_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${SPGLIB_LDFLAGS}"
export CP_LIBS="${SPGLIB_LIBS} \${CP_LIBS}"
EOF
        cat "${BUILDDIR}/setup_spglib" >> $SETUPFILE
fi
cd "${ROOTDIR}"
