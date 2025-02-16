
#include "Hermetis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <time.h>
#include <unordered_map>    
#include <unordered_set>
#include <cmath>
#include "clock.cpp"
#include "data.cpp"
#include "hash_class.cpp"

using namespace std;

#define log_base 1.1

const uint32_t flow_number = 100000;
const char *path = "webget-all-simplify.dat";

struct Webget_Tuple {
    uint64_t id;
    uint64_t timestamp;
};


// vector<data_t> traces;
map<uint64_t, map<uint8_t, uint32_t>> ground_truth_all;
map<uint64_t, map<uint8_t, uint32_t>> ground_truth_large;
map<uint64_t, uint32_t> ground_truth_size;
map<uint32_t, vector<uint64_t>> sort_ground_truth_size;



#define MEMORY 2
uint32_t total_packet_size = 0;
const uint memory = MEMORY * 1000 *1024/8/2;
int Time_threshold = (int)(SLIDING_WINDOW_SIZE) / 3.0
;
sketch<uint32_t, SLOT_NUM, INTERVAL_NUM, HERMETIS_BUCKET_WIDTH, ROW_NUM, SLIDING_WINDOW_SIZE> hist;
    
int main() {
	Webget_Tuple packet;
    vector<data_t> indat;

    unsigned int latestT = 0;
	std::string PATH = "webget-all-simplify.dat";

	int datalen=0;
	std::vector<uint64_t> num1;
	std::vector<uint64_t> num2;

	for (int i=0;i<2;++i)
	{
		std::ifstream file(PATH.c_str());
		while( ! file.eof() )
		{
			uint64_t temp1;
			uint64_t temp2;
			file>>temp1>>temp2;
			if (temp2==0) continue;

			num1.push_back(temp1);
			//num.push_back(temp2 * 123465);
			num2.push_back(temp2);
			
			datalen+=1;
		}
		file.close();
	}

	std::cout<<datalen<<" "<<num1.size()<<" "<<num2.size()<<"\n";
	vector<Webget_Tuple> dataset;
	Webget_Tuple temp;
	for (int i=0;i<datalen;i++)
	{
		temp.id = num1[i];
		temp.timestamp = num2[i];
		dataset.push_back(temp);
	}
	
	std::cout << memory << std::endl;
	// insert data
	uint16_t max_length = 0, min_length = 0xffff;
	uint32_t packet_cnt = 0,packet_CNT = 0;
	double total_insert_time = 0;
	std::cout <<  "Insert Begin" << std::endl;

    queue<Webget_Tuple> dat;

    int cycle = SLIDING_WINDOW_SIZE; // the length of the sliding window, here we use count based sliding window. 
    std::cout << "window size " << cycle << std::endl;

	unsigned int num = 0;
    int Test_num = 0;

	for (auto it = dataset.begin(); it != dataset.end(); it++) {
		if(num % 10000 == 0){
			printf("here \n");
		}
		packet = *it;
		
		int cur_time_tag = ((int)(num / Time_threshold)) % 4;
		// std::cout << "test check point 1" << std::endl;
		uint8_t bid = (log(it->timestamp - 0)  / log(log_base) + 0.5); 
		// std::cout << "bid == "  << (int)bid << std::endl;
		
		new_data_t tmp_key(it->id, bid);

		auto t_a = std::chrono::high_resolution_clock::now();
		hist.insert((it->id), bid, tmp_key, cur_time_tag, num, 1);
		auto t_b = std::chrono::high_resolution_clock::now();
		total_insert_time += std::chrono::duration_cast<std::chrono::microseconds>(t_b - t_a).count();
		packet_CNT++;


		dat.push(packet);
		// std::cout << "test check point 4.0" << std::endl;
		if(ground_truth_all[it->id].find(bid) == ground_truth_all[it->id].end()){
			ground_truth_all[it->id].insert({bid, 1});
		}else{
			ground_truth_all[it->id][bid]++;
		}
		
		if(ground_truth_size.find(it->id) == ground_truth_size.end()){
			ground_truth_size.insert({it->id, 1});
		}else{
			ground_truth_size[it->id]++;
		}
		
		// sort_ground_truth_size[ground_truth_size[it->id] - 1].erase(it->id);
		for(auto iter = sort_ground_truth_size[ground_truth_size[it->id] - 1].begin();iter != sort_ground_truth_size[ground_truth_size[it->id] - 1].end();++iter){
			if(*iter == it->id){
				sort_ground_truth_size[ground_truth_size[it->id] - 1].erase(iter);
				break;
			}
		}
		sort_ground_truth_size[ground_truth_size[it->id]].push_back(it->id);

		if(ground_truth_size[it->id] >= 100){
			for (auto j = ground_truth_all[it->id].begin(); j != ground_truth_all[it->id].end(); ++j) {
				ground_truth_large[it->id][j->first] = j->second;
			}
		}

		max_length = max(max_length, it->timestamp);
		min_length = min(min_length, it->timestamp);
		total_packet_size += it->timestamp;
		packet_cnt++;
		
		while(dat.size() > cycle){

			total_packet_size -= dat.front().timestamp;
			ground_truth_all[dat.front().id][(uint8_t)(log(dat.front().timestamp - 0)  / log(log_base) + 0.5)]--;
			ground_truth_size[dat.front().id]--;
			if(ground_truth_large.find(dat.front().id) != ground_truth_large.end()){
				if(ground_truth_size[dat.front().id] < 100){
					ground_truth_large.erase(ground_truth_large.find(dat.front().id));
				}else{
					for (auto j = ground_truth_all[dat.front().id].begin(); j != ground_truth_all[dat.front().id].end(); ++j) {
						ground_truth_large[dat.front().id][j->first] = j->second;
					}
				}	
			}
			
			for(auto iter = sort_ground_truth_size[ground_truth_size[dat.front().id] + 1].begin();iter != sort_ground_truth_size[ground_truth_size[dat.front().id] + 1].end();++iter){
				if(*iter == dat.front().id){
					sort_ground_truth_size[ground_truth_size[dat.front().id] + 1].erase(iter);
					break;
				}
			}
			sort_ground_truth_size[ground_truth_size[dat.front().id]].push_back(dat.front().id);
			dat.pop();
			packet_cnt--;
		}
		num++;

		

		if( ((int)packet.timestamp - (int)indat[0].timestamp) % 2*cycle == 0){
		
			// std::cout << "test index begin" << std::endl;
			hist.statistics.pointLarge = 0;
			hist.statistics.pointAll = 0;
			hist.statistics.histogramLarge = 0;
			hist.statistics.histogramAll = 0;
			hist.statistics.pointLargeARE = 0;
			hist.statistics.pointAllARE = 0;
			hist.statistics. histogramLargeARE = 0;
			hist.statistics.histogramAllARE = 0;
			hist.statistics.cardinalityRE = 0;
			hist.statistics.entropyRE = 0;
			hist.statistics.pointChange = 0;
			hist.statistics.histogramChange = 0;

            latestT = packet.timestamp;
			double are = 0;
			uint32_t bucket_large = 0;
			packet_CNT = 0;
			total_insert_time = 0;
			
			for (auto it = ground_truth_large.begin(); it != ground_truth_large.end(); ++it) {
				if(hist.heavyPart.isExist((Key_t)it->first.str)){
					for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {			// true histogram
						uint32_t result = hist.pointQuery((Key_t)it->first.str, it2->first, num);
						
						are = (it2->second != 0) ? (abs((double)result - it2->second) / it2->second) : 0;
						hist.statistics.pointLargeARE += are;
						if (are < POINT_ARETHRESHOLD) {
							hist.statistics.pointLarge++;
						}
						bucket_large++;
					}
				}
				
			}

			// point query of all flows
			uint32_t bucket_all = 0;
			int complex_flow_count = 0;
			for (auto it = ground_truth_all.begin(); it != ground_truth_all.end(); ++it) {
				if(it->second.size()>5) complex_flow_count++;
				for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {			// true histogram
					uint32_t result = hist.pointQuery((Key_t)it->first.str, it2->first, num);
					hist.statistics.pointAllARE += are;
					are = (it2->second != 0) ? (abs((double)result - it2->second) / it2->second) : 0;
					if (are < POINT_ARETHRESHOLD) {
						hist.statistics.pointAll++;
					}
					
					bucket_all++;
				}
			}
			std::cout << "complex_flow_count is " << complex_flow_count << std::endl;

			//histogram query of large flows
			for (auto it = ground_truth_large.begin(); it != ground_truth_large.end(); ++it) {
				double totalerr = 0;
				uint32_t totalcnt = 0;
				Hist_t result = hist.histogramQuery((Key_t)it->first.str,(unsigned long long int)num);
				for (int i = 0; i  < BUCKET_NUM; ++i) {
					if (it->second.find(i) != it->second.end()) {
						totalerr += fabs((double)result.at(i) - it->second[i]);
						totalcnt += it->second[i];
					}
					else {
						totalerr += fabs(result.at(i));
					}
				}
			are = (totalcnt != 0) ? (totalerr / totalcnt) : 0;
				hist.statistics.histogramLargeARE += are;
				if (are < HISTOGRAM_ARETHRESHOLD) {
					hist.statistics.histogramLarge++;
				}
			}

			// histogram query of all flows
			for (auto it = ground_truth_all.begin(); it != ground_truth_all.end(); ++it) {
				double totalerr = 0;
				uint32_t totalcnt = 0;
				Hist_t result = hist.histogramQuery((Key_t)it->first.str, (unsigned long long int)num);
				for (int i = 0; i  < BUCKET_NUM; ++i) {
					if (it->second.find(i) != it->second.end() && it->second[i]) {
						totalerr += fabs((double)result.at(i) - it->second[i]);
						totalcnt += it->second[i];
					}
					else {
						totalerr += fabs(result.at(i));
					}
				}
				are = (totalcnt != 0) ? (totalerr / totalcnt) : 0;
				// std::cout << "totalerr == " << totalerr << ", totalcnt == " << totalcnt << ", are == "  << are << std::endl;
				hist.statistics.histogramAllARE += are;
				if (are < HISTOGRAM_ARETHRESHOLD) {
					hist.statistics.histogramAll++;
				}
			
			}

			// cardinality
			std::map<int, int> ground_truth_distribution;
			for (auto it = ground_truth_all.begin(); it != ground_truth_all.end(); it++) {
				const int value = it->second.size();
				// std::cout <<"value is " << value << std::endl;
				if (ground_truth_distribution.find(value) == ground_truth_distribution.end())
					ground_truth_distribution.emplace(value, 1);
				else
					ground_truth_distribution[value] += 1;
			}  

			uint32_t cardinality_truth = get_cardinality(ground_truth_all);
			uint32_t cardinality_estimate = hist.get_cardinality(cur_time_tag);
			hist.statistics.cardinalityRE = fabs((double)cardinality_truth - cardinality_estimate) / cardinality_truth;
			std::cout << "Cardinality(estimate, truth):" << cardinality_estimate << " " << cardinality_truth << std::endl;

			// entropy
			double entropy_truth = get_entropy(ground_truth_all);
			double entropy_estimate = hist.get_entropy(cur_time_tag);
			hist.statistics.entropyRE = fabs((double)entropy_truth - entropy_estimate) / entropy_truth;
   			
			std::cout << "Point query for large flows (proportion of buckets): " <<
				(double)hist.statistics.pointLarge / bucket_large << " " << hist.statistics.pointLargeARE / bucket_large << std::endl;
			std::cout << "Point query for all flows (proportion of buckets): " <<
				(double)hist.statistics.pointAll / bucket_all << " " << hist.statistics.pointAllARE / bucket_all << std::endl;
			std::cout << "Histogram query for large flows (proportion of flows): " <<
				(double)hist.statistics.histogramLarge / ground_truth_large.size() << " " << hist.statistics.histogramLargeARE / ground_truth_large.size() << std::endl;
			std::cout << "Histogram query for all flows (proportion of flows): " <<
				(double)hist.statistics.histogramAll / ground_truth_all.size() << " " << hist.statistics.histogramAllARE / ground_truth_all.size() << std::endl;
			std::cout << "Cardinality relative error: " << hist.statistics.cardinalityRE << std::endl;
			std::cout << "Entropy relative error: " << hist.statistics.entropyRE << std::endl;
			std::cout << "Bandwidth: " << (double)hist.bandwidth / total_packet_size << " " << (double)CM_DEPTH_HIST * CM_WIDTH_HIST * sizeof(uint16_t) / total_packet_size << " " << (double)(hist.bandwidth + CM_DEPTH_HIST * CM_WIDTH_HIST * sizeof(uint16_t)) / total_packet_size << std::endl;
			std::cout << "nonempty : " << hist.heavyPart.nonempty << std::endl;
			cout << "Memory usage: " << (double)hist.get_memory_usage() / 1024 / 1024 << " MB" << endl;
	}
		
}
	return 0;
}


