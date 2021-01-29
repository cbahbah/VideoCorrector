cmake_minimum_required (VERSION 2.6)

project (VideoEditorCB)


find_package( OpenCV REQUIRED )

set(VideoEditorCB_SOURCES
	main.cpp
    VideoSequence.cpp
)

set(VideoEditorCB_HEADERS
	VideoSequence.h
)


add_executable(VideoEditorCB ${VideoEditorCB_SOURCES} ${VideoEditorCB_HEADERS})

target_link_libraries( VideoEditorCB
    ${OpenCV_LIBS})