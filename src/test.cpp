
#include <iostream>
#include <vector>
#include <algorithm>


static std::vector<std::string> extractDigits(const std::string& inputString);

int countElementOccurrences(const std::vector<int>& vec, int targetElement) {
    int occurrences = std::count(vec.begin(), vec.end(), targetElement);
    return occurrences;
}

int main(int argc, char* argv[]){
	std::vector<int> numbers = {1, 2, 3, 4, 2, 5, 2, 6, 2, 7};

    int targetElement = 2;

    int occurrences = countElementOccurrences(numbers, targetElement);

    std::cout << "Element " << targetElement << " occurs " << occurrences << " times." << std::endl;

    return 0;

}

std::vector<std::string> extractDigits(const std::string& inputString) {
    std::vector<std::string> strs;
    std::string tempString = "";

    for (char c : inputString) {
        if (c != '_') {
            tempString += c;
        } else {
            if (!tempString.empty()) {
                strs.push_back(tempString);
                tempString.clear();
            }
        }
    }
    if (!tempString.empty()) {
        strs.push_back(tempString);
    }

    return strs;
}
