#include <TApplication.h>
#include "FaserCalDisplay.h"

int main(int argc, char **argv) {
    TApplication app("FaserCalDisplayApp", &argc, argv);
    display::FaserCalDisplay display;
    display.GetEventDisplay();
    //display.GetDetector();
    app.Run();
    return 0;
}

