#include <iostream>

#include "asvh.h"

using namespace std;

int main()
{
    cout << "The Recesions Analyses" << endl;

    asvh recesions("./settings/settings_file");
    recesions.get_recessions_data();
    recesions.initialize_all_POpulation();
    recesions.run_ensemble();


    return 0;
}
