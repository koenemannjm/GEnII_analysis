#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include <TString.h>
#include <TMap.h>
#include <TObjString.h>


TMap* gConfigMap = nullptr;
std::string CONFIG_PATH = "/work/halla/sbs/koeneman/GEnII/GEnII_analysis/config/";

void readConfig(const std::string& config_filename) {

  std::string config_filepath = CONFIG_PATH + config_filename;

  if (!gConfigMap) {
    gConfigMap = new TMap();
    gConfigMap->SetOwnerKeyValue(kTRUE, kTRUE);
  } 

  std::ifstream configFile(config_filepath);
  if (!configFile.is_open()) {
    std::cerr << "Error >> Config file does not exist: " << config_filepath << std::endl;
    return;
  }

  std::string line;
  int lineNumber = 0;
  int parsedCount = 0;

  while (std::getline(configFile, line)) {
    lineNumber++;

    if (line.empty() || line[0] == '#') continue;

    size_t start = line.find_first_not_of(" \t");
    size_t end = line.find_last_not_of(" \t");
    if (start == std::string::npos) continue;
    line = line.substr(start, end - start + 1);

    size_t spacePos = line.find(' ');
    if (spacePos == std::string::npos) {
      std::cerr << "Warning >> Line " << lineNumber
		<< " has no space seperator: " << line << std::endl;
      continue;
    }

    std::string variableStr = line.substr(0, spacePos);
    std::string valueStr = line.substr(spacePos+1);
    valueStr.erase(valueStr.find_last_not_of(" \t") + 1);

    gConfigMap->Add(new TObjString(variableStr.c_str()),
		    new TObjString(valueStr.c_str()));
    parsedCount++;
    
  }
  configFile.close();
  
}

TString getConfigString(const char* key, const char* defaultValue = "") {
  if (!gConfigMap) {
    return TString(defaultValue);    
  }

  TObjString* keyObj = new TObjString(key);
  TObjString* valueObj = (TObjString*)gConfigMap->GetValue(keyObj);
  delete keyObj;

  if (valueObj) {
    return valueObj->GetString();
  }

  return TString(defaultValue);
}

int getConfigInt(const char* key, int defaultValue = 0) {
  TString value = getConfigString(key);
  if (value.IsNull()) return defaultValue;
  return value.Atoi();
}

double getConfigDouble(const char* key, double defaultValue = 0.0) {
  TString value = getConfigString(key);
  if (value.IsNull()) return defaultValue;
  return value.Atof();
}

void printConfig() {
  if (!gConfigMap) {
    std::cout << "Problem! No config loaded!" << std::endl;
    return;
  }

  std::cout << "\n=== Config File ===" << std::endl;
  TIterator* iter = gConfigMap->MakeIterator();
  TObjString* key;
  int lineNumber = 0;
  while ((key = (TObjString*)iter->Next())) {
    lineNumber++;
    TObjString* value = (TObjString*)gConfigMap->GetValue(key);
    std::cout << " [" << lineNumber << "] "
	      << key->GetString() << " = "
	      << value->GetString() << std::endl;
  }
  delete iter;
}

void clearConfig() {
  if (gConfigMap) {
    delete gConfigMap;
    gConfigMap = nullptr;
  }
}

// === TESTING CONFIG PARSER ===

void testConfig(const std::string& config_filename) {

  readConfig(config_filename);

  printConfig();
}
