#!/bin/sh

# made using https://github.com/JulianoGianlupi/cc3d-nanoHub-tool-setuper
#made with help from Steven Clark

current_directory=$(pwd)

echo ${current_directory}

PARENT_DIR="$(dirname "$(pwd)")"

echo ${PARENT_DIR}

scriptPath=$0

binDir=$(dirname ${scriptPath})

export NANOHUB_SIM=${binDir}

TOOL_HOME=$(dirname ${binDir})

echo ${TOOL_HOME}

# necessary to enforce standard convention for numeric values specification on non-English OS
export LC_NUMERIC="C.UTF-8"


CC3D_TOOL_NAME=cctrinity
CC3D_TOOL_FILENAME=cctrinity.cc3d
# CC3D_DEMO_PREFIX=Demos/BookChapterDemos_ComputationalMethodsInCellBiology/VascularTumor/
USER_DEMO_PATH=$HOME/CompuCell3D_Demos/${CC3D_TOOL_NAME}
echo ${USER_DEMO_PATH}
MAIN_CC3D_FILES=${TOOL_HOME}/main/.
echo ${MAIN_CC3D_FILES}


mkdir -p ${USER_DEMO_PATH}

mkdir -p ${USER_DEMO_PATH}/output/

cp -rf ${MAIN_CC3D_FILES} ${USER_DEMO_PATH}/ #not checking for existance incase there's been an update

${COMPUCELL3D_HOME}/compucell3d.sh -i ${USER_DEMO_PATH}/${CC3D_TOOL_FILENAME} -o ${USER_DEMO_PATH}/output/ -f 10001 > /dev/null 2>&1 
