/*
 * VideoSequence.cpp
 *
 *  Created on: 26 janv. 2021
 *      Author: chahrazadebahbah
 */

#include <iostream>
#include <fstream>

#include "VideoSequence.h"
#include "VideoSequence.h"
#include <numeric>
#include <stdlib.h>

#include <math.h>


#include <opencv2/opencv.hpp>
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"

using namespace cv;
using namespace cv::xfeatures2d;


namespace
{
typedef std::vector<std::pair<int, double>> IndexesMeanSSIMs;
typedef std::pair<int, double> FrameSSIM;
typedef std::vector<FrameSSIM> AllFramesSSIM;
typedef std::pair<int, AllFramesSSIM>  Ref_AllFramesSSIM;
std::vector<Ref_AllFramesSSIM> framesSSIM;


Scalar computeSSIM_index(const Mat &i1, const Mat &i2)
{
	//Taken from https://docs.opencv.org/master/dd/d3d/tutorial_gpu_basics_similarity.html
	const double C1 = 6.5025, C2 = 58.5225;

	/***************************** INITS **********************************/
	int d = CV_32F;

	Mat I1, I2;

	i1.convertTo(I1, d);           // cannot calculate on one byte large values
	i2.convertTo(I2, d);
	Mat I2_2 = I2.mul(I2);        // I2^2
	Mat I1_2 = I1.mul(I1);        // I1^2
	Mat I1_I2 = I1.mul(I2);        // I1 * I2

	/*************************** END INITS **********************************/

	Mat mu1, mu2;   // PRELIMINARY COMPUTING
	GaussianBlur(I1, mu1, Size(11, 11), 1.5);
	GaussianBlur(I2, mu2, Size(11, 11), 1.5);

	Mat mu1_2 = mu1.mul(mu1);
	Mat mu2_2 = mu2.mul(mu2);
	Mat mu1_mu2 = mu1.mul(mu2);


	Mat sigma1_2, sigma2_2, sigma12;

	GaussianBlur(I1_2, sigma1_2, Size(11, 11), 1.5);
	sigma1_2 -= mu1_2;
	GaussianBlur(I2_2, sigma2_2, Size(11, 11), 1.5);
	sigma2_2 -= mu2_2;
	GaussianBlur(I1_I2, sigma12, Size(11, 11), 1.5);
	sigma12 -= mu1_mu2;

	Mat t1, t2, t3;
	t1 = 2 * mu1_mu2 + C1;
	t2 = 2 * sigma12 + C2;
	t3 = t1.mul(t2);     // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

	t1 = mu1_2 + mu2_2 + C1;
	t2 = sigma1_2 + sigma2_2 + C2;
	t1 = t1.mul(t2);   // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

	Mat ssim_map;
	divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;
	Scalar mssim = mean(ssim_map); // mssim = average of ssim map

	return mssim;
}

void storeSSIM(std:: vector<std::pair<int,cv::Mat>>& allFrames, std::string file1, std::string file2)
{
	//Store in text file (heavy computation)
	ofstream storedSSIM(file1);
	ofstream storedMeanSSIM(file2);

	Scalar SSIM;
	std::vector<std::pair<int,double>> framesSSIM;



	for (auto& af : allFrames)
	{
		storedSSIM << af.first<< " ";
		storedMeanSSIM << af.first<< " ";

		double finalSSIM = 0.0;
		for (auto& afBis : allFrames)
		{
			if (af.first == afBis.first)
			{
				SSIM = 0.0;
			}
			else
			{
				SSIM = computeSSIM_index(af.second, afBis.second);
			}
			auto meanSSIM = cv::sum(SSIM)/3.0;
			finalSSIM += meanSSIM[0];
			storedSSIM << afBis.first<< " "<<meanSSIM[0]<<" ";


		}
		storedSSIM << " \n";
		storedMeanSSIM << finalSSIM/allFrames.size()<<" \n";


		framesSSIM.push_back(std::make_pair(af.first, finalSSIM/allFrames.size()));



	}
	for (auto& p : framesSSIM)
	{
		std::cout<<"Frames " <<p.first<<"  MeanSSIM "<<p.second<<std::endl;

	}
}

void storeSURF(std:: vector<std::pair<int,cv::Mat>>& allFrames)
{

    //-- Step 1: Detect the keypoints using SURF Detector
	int minHessian = 400;
	Ptr<SURF> detector = SURF::create( minHessian );
	typedef std::vector<KeyPoint> KeyPoints;
	KeyPoints keypoints;
	std::vector<std::pair<int,Mat>> allDescriptors;


	Mat descriptors;
	for (auto& aF : allFrames)
	{
		detector->detectAndCompute( aF.second, noArray(), keypoints, descriptors );
		allDescriptors.push_back(std::make_pair(aF.first, descriptors));
	}

    //-- Step 2: Matching descriptor vectors with a FLANN based matcher
    // Since SURF is a floating-point descriptor NORM_L2 is used
    Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(DescriptorMatcher::FLANNBASED);
    typedef std::vector< std::vector<DMatch> > KnnMatches;
    KnnMatches knn_matches;
    std::vector<KnnMatches> allKnn_matches;
	for (auto& aPair : allDescriptors)
	{
		for (auto& aPairBis : allDescriptors)
		{
		    matcher->knnMatch( aPair.second, aPairBis.second, knn_matches, 2 );
		    allKnn_matches.push_back(knn_matches);

		}
	}


	//-- Filter matches using the Lowe's ratio test
	const float ratio_thresh = 0.7f;
	std::vector<std::vector<DMatch>> allGood_matches;

	for (auto& aMatch : allKnn_matches)
	{
		std::vector<DMatch> good_matches;
		for (size_t i = 0; i < aMatch.size(); ++i)
		{
			if (aMatch[i][0].distance < ratio_thresh * aMatch[i][1].distance)
			{
				good_matches.push_back(aMatch[i][0]);
				std::cout<<"good matches "<<aMatch[i][0].distance<<std::endl;

			}
			allGood_matches.push_back(good_matches);


		}

	}


}

IndexesMeanSSIMs readMeanSSIM()
{
	string line;
	ifstream myfile("meanSSIM.txt");
	std::vector<double> meanSSIMs;
	IndexesMeanSSIMs indexesMeanSSIMs;

	std::string frameIndex, meanSSIM;
	int count = 0;
	if (myfile.is_open())
	{
	    while ( getline (myfile,line) )
	    {
	    	std::stringstream ss(line);
	    	std::getline(ss,frameIndex,' ');
	    	std::getline(ss,meanSSIM,' ');
	    	int iFrameIndex = stoi(frameIndex);
	    	double dMeanSSIM = stod(meanSSIM);
	    	indexesMeanSSIMs.push_back(std::make_pair(iFrameIndex, dMeanSSIM));
	    	meanSSIMs.push_back(dMeanSSIM);
	    	++count;
	    }
	    myfile.close();
	  }

	  else
	  {
		  cout << "Unable to open file";
	  }


	std::sort(meanSSIMs.begin(),meanSSIMs.end());
	auto sum = std::accumulate(meanSSIMs.begin(),meanSSIMs.end(),0.0);
	double epsMoy = sum/count;

	auto sumVar = std::accumulate(meanSSIMs.begin(),meanSSIMs.end(),0.0,
			[&epsMoy](double moy, double aM2)
			{
				moy += (aM2 - epsMoy)*(aM2 - epsMoy);
				return moy;
			});

	double epsVar = sumVar/(count-1);
	double eps = epsVar + epsMoy;

	int size = meanSSIMs.size();
	int medianIndex = (size + 1)/2;
	double epsMedian = meanSSIMs[ medianIndex];
	IndexesMeanSSIMs selectedIndexesMeanSSIMs;
	//more tests should be performed to define a proper threshold
	for (auto& aPair : indexesMeanSSIMs)
	{
		if (aPair.second > epsMedian)
		{
			std::cerr<<"frameOK "<<aPair.first<< " "<<aPair.second<<std::endl;
			selectedIndexesMeanSSIMs.push_back(aPair);
		}
		else
		{
			std::cerr<<"frameKO "<<aPair.first<< " "<<aPair.second<<std::endl;

		}


	}

	return selectedIndexesMeanSSIMs;

}

std::vector<Ref_AllFramesSSIM> readSSIM(std::string file1)
{
	string line;
	ifstream myfile(file1);
	std::vector<Ref_AllFramesSSIM> framesSSIM;


	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{

			std::stringstream ss(line);
			std::string tmp;
			std::vector<std::string> tokens;
			while (getline(ss, tmp, ' '))
			{
				tokens.push_back(tmp);
			}

			int initFrame = stoi(tokens[0]);
			AllFramesSSIM tmpFrameSSIM;
			int pos = 0;
			int frame = 0;
			double ssim = 0.0;

			std::vector<double> tmpValues;
			std::vector<int> tmpFrs;
			for (int i = 1; i< tokens.size()-1; ++i)
			{

				if (i%2 != 0)
				{
					frame = stoi(tokens[i]);
					tmpFrs.push_back(frame);
					++pos;
				}

				else
				{

					ssim = stod(tokens[i]);
					tmpValues.push_back(ssim);
					++pos;
				}
				if (pos == 2)
				{

					tmpFrameSSIM.push_back(std::make_pair(tmpFrs.front(), tmpValues.front()));
					tmpFrs.erase(tmpFrs.begin());
					tmpValues.erase(tmpValues.begin());

					pos = 0;

				}

			}


			framesSSIM.push_back(std::make_pair(initFrame,tmpFrameSSIM));


		}
		myfile.close();

	}




	return framesSSIM;

}

std::vector<int> orderSequence(std::vector<int> indexesGoodFrames)
{
	auto res = readSSIM("SSIM_filtered.txt");

	int refIndex = res.front().first;

	std::vector<int> sortedFrames;
	sortedFrames.push_back(refIndex);
	std::vector<double> valuesGoodSSIM;
	indexesGoodFrames.erase(indexesGoodFrames.begin());



	while (!indexesGoodFrames.empty())
	{


		auto find = std::find_if(res.begin(), res.end(),
				[&](Ref_AllFramesSSIM& aRes)
				{
			return aRes.first == refIndex;
				});

		if (find != res.end())
		{
			//res.erase(find);

			std::vector<std::pair<int,double>> aVec;
			for (auto& v : find->second)
			{
				auto checkFind = std::any_of(sortedFrames.begin(), sortedFrames.end(),
						[&](int aFrame)
						{
					return v.first == aFrame;
						});


				if (!checkFind)
				{
					aVec.push_back(v);
				}

			}


			auto max = std::max_element(aVec.begin(), aVec.end(),
					[&](std::pair<int, double>& aVal1, std::pair<int, double>& aVal2)
					{
						return aVal1.second < aVal2.second;
					});
			sortedFrames.push_back(max->first);

			valuesGoodSSIM.push_back(max->second);

			refIndex = max->first;

			std::vector<int>::iterator indexToDelete = std::find(indexesGoodFrames.begin(),
					indexesGoodFrames.end(), max->first);

			res.erase(find);

			indexesGoodFrames.erase(indexToDelete);





		}


	}
	return sortedFrames;

}


}


VideoSequence::VideoSequence()
{
}


VideoSequence::~VideoSequence() {
	// TODO Auto-generated destructor stub
}





void VideoSequence::loadSequence(std::string aInputVideo)
{
	cv::VideoCapture videoCap(aInputVideo);
	if (videoCap.isOpened())
	{
		std::cerr<<"Open Video "<<aInputVideo<<std::endl;



		std:: vector<std::pair<int,cv::Mat>> allFrames;
		int indexFrame = 0;
		while (videoCap.isOpened())
		{
			Mat frame;
			videoCap>>frame;

			if (frame.empty())
			    break;

			allFrames.push_back(std::make_pair(indexFrame, frame));
			++indexFrame;
			//cv::imshow("Frame", frame);
		}

		videoCap.release();

		//Compute and store meanSSIM in file text
		//storeSSIM(allFrames, "SSIM.txt", "meanSSIM.txt");


		//select good frames (filter outliers)
		IndexesMeanSSIMs selectedframes = readMeanSSIM();
		std:: vector<std::pair<int,cv::Mat>> goodFrames;
		std::vector<int> indexesGoodFrames;

		for (auto& aPair : allFrames )
		{
			auto found = std::any_of(selectedframes.begin(), selectedframes.end(),
					[&](const std::pair<int, double>& aP)
					{
						return aPair.first == aP.first;
					});
			if (found)
			{
				goodFrames.push_back(std::make_pair(aPair.first, aPair.second));
				indexesGoodFrames.push_back(aPair.first);
			}


		}

		//Compute SSURF between selected frames
		//storeSURF(goodFrames);

		double frame_fps = videoCap.get(CAP_PROP_FPS);;
		int frame_width = videoCap.get(CAP_PROP_FRAME_WIDTH);
		int frame_height = videoCap.get(CAP_PROP_FRAME_HEIGHT);


		cv::VideoWriter video("/Users/chahrazadebahbah/eclipse-workspace/VideoEditorCB/filtered.mp4",
				VideoWriter::fourcc('m', 'p', '4', 'v'),frame_fps
					, Size(frame_width,frame_height));

		for (auto& aP : allFrames)
		{
			video.write(aP.second);
		}

		video.release();




		//Compute SSIM between filtered frames
		//storeSSIM(goodFrames, "SSIM_filtered.txt", "meanSSIM_filtered.txt");



		//Order sequence, for now we start from the first frame
		std::vector<int> finalFrames = orderSequence(indexesGoodFrames);
		std::vector<std::pair<int,Mat>> finalOragnisedFrames;
		std::transform(finalFrames.begin(), finalFrames.end(),
				 std::back_inserter(finalOragnisedFrames),
				[&](int aIndex)
				{
					auto found = std::find_if(goodFrames.begin(), goodFrames.end(),
							[&](const std::pair<int, Mat>& aP)
							{
								return aP.first == aIndex;
							});

					if (found != goodFrames.end())
					{
						return std::make_pair(aIndex, found->second);
					}
				});




		//Write final video

		cv::VideoWriter videoFinal("/Users/chahrazadebahbah/eclipse-workspace/VideoEditorCB/finalVideo.mp4",
				VideoWriter::fourcc('m', 'p', '4', 'v'),frame_fps, Size(frame_width,frame_height));

		for (auto& aP : allFrames)
		{
			videoFinal.write(aP.second);
		}
		videoFinal.release();

		cv::destroyAllWindows();




	}
	else
	{
		std::cerr<<"Close video to begin "<<std::endl;
	}



}
