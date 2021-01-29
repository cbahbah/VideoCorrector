/*
 * VideoSequence.h
 *
 *  Created on: 26 janv. 2021
 *      Author: chahrazadebahbah
 *      "Class to load video"
 */

#ifndef VIDEOSEQUENCE_H_
#define VIDEOSEQUENCE_H_

#include <string>

using namespace std;


class VideoSequence {
public :
	VideoSequence();
	virtual ~VideoSequence();

	void loadSequence(std::string inputVideo);
	void computeSSIM();


private :
};




#endif /* VIDEOSEQUENCE_H_ */
