#设置cmake版本号
cmake_minimum_required(VERSION 3.4.0)

##########################################################################
#设置工程名称
project(MonteCarlo)

#设置工程包含当前目录，非必须
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

##########################################################################
#配置CUDA
# find_package(CUDA REQUIRED)

##########################################################################
#配置VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

##########################################################################
#配置Qt
#设置自动生成moc文件，一定要设置
set(CMAKE_AUTOMOC ON) 

#查找需要的Qt库文件，最好每一个库都要写，Qt也会根据依赖关系自动添加
find_package(Qt5Widgets)  
find_package(Qt5Core)
find_package(Qt5Gui)
find_package(Qt5LinguistTools)

##########################################################################
#配置项目源码
#查找当前文件夹中的所有源代码文件，也可以通过Set命令将所有文件设置为一个变量
FILE(GLOB SRC_FILES "./*.cpp") 

#同理，查找当前文件夹中所有的头文件
FILE(GLOB HEAD_FILES "./*.h") 

#同理，查找当前文件夹中所有的ui文件
FILE(GLOB UI_FILES "./UI/*.ui") 
#同理，查找当前文件夹中所有的qrc文件
FILE(GLOB RCC_FILES "./UI/*.qrc") 
#同理，查找当前文件夹中所有的qss文件
FILE(GLOB QSS_FILES "./UI/qss/*.qss") 
#同理，查找当前文件夹中所有的ts文件
FILE(GLOB TS_FILES "./UI/language/*.ts") 

FILE(GLOB MONTE_CARLO_FILES "./MonteCarlo/*.cpp" "./MonteCarlo/*.h")

#通过Ui文件生成对应的头文件，一定要添加
qt5_wrap_ui(WRAP_FILES ${UI_FILES})
qt5_add_resources(QRC_FILES ${RCC_FILES})

##########################################################################
#源码分组
#将ui文件和生成文件整理在一个文件夹中，非必须
source_group("Ui" FILES ${UI_FILES})
source_group("Ui_qrc" FILES ${RCC_FILES})
source_group("Ui_qss" FILES ${QSS_FILES})
source_group("Ui_language" FILES ${TS_FILES})

source_group("Source" FILES ${SRC_FILES})
source_group("Header" FILES ${HEAD_FILES})

source_group("MonteCarlo" FILES ${MONTE_CARLO_FILES})
##########################################################################
#创建工程

#翻译文件
qt5_add_translation(QM_FILES ${TS_FILES})

#lupdate翻译文件
add_custom_target (lupdate_ts
COMMAND lupdate -recursive ${CMAKE_CURRENT_SOURCE_DIR} -ts ${CMAKE_CURRENT_SOURCE_DIR}/UI/language/language_zh.ts
)
#lrelease翻译文件
add_custom_target (lrelease_ts
COMMAND lrelease ${CMAKE_CURRENT_SOURCE_DIR}/UI/language/language_zh.ts -qm ${CMAKE_CURRENT_SOURCE_DIR}/UI/language/language_zh.qm
)
#让lrelease依赖于lupdate，生成lrelease是先生成lupdate
#add_dependencies(lrelease_ts lupdate_ts)

add_executable(${PROJECT_NAME} ${SRC_FILES} ${HEAD_FILES} ${QRC_FILES} ${WRAP_FILES} ${QSS_FILES} ${TS_FILES} ${MONTE_CARLO_FILES})

#让本项目依赖于lrelease，生成时先生成项目lrelease
#add_dependencies(${PROJECT_NAME} lrelease_ts)

##########################################################################
#添加Qt5依赖项
target_link_libraries(${PROJECT_NAME} Qt5::Widgets Qt5::Core Qt5::Gui)
##########################################################################
#添加VTK依赖项
if(VTK_LIBRARIES)
  target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
else()
  target_link_libraries(${PROJECT_NAME} vtkHybrid vtkWidgets)
endif()

#添加VTK依赖项
if(SVR_LIBRARY)
  target_link_libraries(${PROJECT_NAME} ${SVR_LIBRARIES})
endif()

