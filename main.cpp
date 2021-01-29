#include <iostream>
#include "VideoSequence.h"

int main(int argc, char **argv) {
	std::cout<<"Fix corrupted video"<<std::endl;
	std::string pathVideo = "/Users/chahrazadebahbah/eclipse-workspace/VideoEditorCB/corrupted_video.mp4";

	VideoSequence videoSeq;
	videoSeq.loadSequence(pathVideo);
	return 0;
}


