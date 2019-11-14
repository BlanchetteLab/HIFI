HIFI: src/HIFI.cpp src/HIFI_options.cpp
	g++ -O4 -o src/HIFI src/HIFI.cpp src/HIFI_options.cpp -lm
callPeaks: src/callPeaks.cpp
	g++ -O4 -o src/callPeaks src/callPeaks.cpp -lm
