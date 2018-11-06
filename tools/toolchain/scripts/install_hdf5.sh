#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

hdf5_ver=${hdf5_ver:-1.10.4}

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_hdf5=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_hdf5" ] && rm "${BUILDDIR}/setup_hdf5"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_hdf5" in
    __INSTALL__)
        echo "==================== Installing hdf5 ===================="
        pkg_install_dir="${INSTALLDIR}/hdf5-${hdf5_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "hdf5-${hdf5_ver} is already installed, skipping it."
        else
            if [ -f hdf5-${hdf5_ver}.tar.bz2 ] ; then
                echo "hdf5-${hdf5_ver}.tar.bz2 is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/hdf5-${hdf5_ver}.tar.bz2
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d hdf5-${hdf5_ver} ] && rm -rf hdf5-${hdf5_ver}
            tar xf hdf5-${hdf5_ver}.tar.bz2
            cd hdf5-${hdf5_ver}
            ./configure \
                --prefix="${pkg_install_dir}" \
                --libdir="${pkg_install_dir}/lib" \
                --enable-shared \
                > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi

        HDF5_CFLAGS="-I${pkg_install_dir}/include"
        HDF5_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding hdf5 from system paths ===================="
        check_command pkg-config --modversion hdf5
        add_include_from_paths HDF5_CFLAGS "hdf5.h" $INCLUDE_PATHS
        add_lib_from_paths HDF5_LDFLAGS "libhdf5.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        ;;
esac
if [ "$with_hdf5" != "__DONTUSE__" ] ; then
    HDF5_LIBS="-lhdf5"
    if [ "$with_hdf5" != "__SYSTEM__" ] ; then
    cat <<EOF > "${BUILDDIR}/setup_hdf5"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    fi
    cat <<EOF > "${BUILDDIR}/setup_hdf5"
export HDF5_CFLAGS="${HDF5_CFLAGS}"
export HDF5_LDFLAGS="${HDF5_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(IF_OMP(-D__HDF5|)|)"
export CP_CFLAGS="\${CP_CFLAGS} ${HDF5_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${HDF5_LDFLAGS}"
####################################################
#
# inlcude hdf5 only if sirius is activated and build
# depends them on mpi and omp
#
####################################################

export CP_LIBS="IF_MPI(IF_OMP(${HDF5_LIBS}|)|) \${CP_LIBS}"
EOF
    cat "${BUILDDIR}/setup_hdf5" >> $SETUPFILE
fi
cd "${ROOTDIR}"
