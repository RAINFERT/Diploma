#pragma once

#include "ComponentDatabase.h"

#include <string>

class ComponentDatabaseLoader {
public:
    static ComponentDatabase loadFromCsv(
        const std::string& path
        );
};
