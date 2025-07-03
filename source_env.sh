DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
export PYTHONPATH=${PYTHONPATH}:${DIR}
export PATH=${PATH}:${DIR}/bin
export PHENIX_TRUST_OTHER_ENV=1
