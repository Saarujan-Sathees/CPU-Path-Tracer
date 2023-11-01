#ifndef JSON_STRUCTURE
#define JSON_STRUCTURE

#include <vector>
#include <algorithm>
#include "Vector.h"

struct Types;

class JSON {
    private:
    std::vector<std::string> keys;
    std::vector<Types*> values;

    public:
    JSON() {
        keys = std::vector<std::string>(0);
        values = std::vector<Types*>(0);
    }

    void addPair(const std::string& key, Types* value) {
        keys.push_back(key);
        values.push_back(value);
    }

    bool contains(const std::string& key) const {
        std::string lowercaseKey = "";
        for (char c : key) {
            lowercaseKey += tolower(c);
        }

        return std::find(keys.begin(), keys.end(), lowercaseKey) != keys.end();
    }

    Types* get(const std::string& key) const {
        std::string lowercaseKey = "";
        for (char c : key) {
            lowercaseKey += tolower(c);
        }

        int index = std::find(keys.begin(), keys.end(), lowercaseKey) - keys.begin();
        if (index == keys.end() - keys.begin()) 
            return nullptr;
        
        return values[index];
    }

    const std::vector<std::string> getKeys() const {
        return keys;
    }

    const std::vector<Types*> getValues() const {
        return values;
    }
};

struct Types {
    double dblValue = INFINITY;
    Vector* vecValue = nullptr; 
    JSON jsonValue;
    bool boolValue = false;
    std::string strValue = "";
    std::vector<Types*> arrValue = std::vector<Types*>(0);

    Types() {}
    Types(double dbl) { dblValue = dbl; }
    Types(Vector* vec) { vecValue = vec; }
    Types(JSON json) { jsonValue = json; }
    Types(bool boolean) { boolValue = boolean; }
    Types(std::string str) { strValue = str; }
    Types(std::vector<Types*> arr) { arrValue = arr; }
};


#endif