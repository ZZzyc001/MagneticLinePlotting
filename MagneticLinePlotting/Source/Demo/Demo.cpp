#include "ArgsParser.h"
#include "MagneticSandBuilder.h"
#include "Simulator.h"

#include <iostream>
using namespace PhysX;

inline std::unique_ptr<ArgsParser> BuildArgsParser() {
    auto parser = std::make_unique<ArgsParser>();
    parser->addArgument<std::string>("output", 'o', "the output directory", "output");
    parser->addArgument<int>("dim", 'd', "the dimension of the simulation", 3);
    parser->addArgument<int>("test", 't', "the test case index", 12);
    parser->addArgument<uint>("begin", 'b', "the begin frame (including)", 0);
    parser->addArgument<uint>("end", 'e', "the end frame (excluding)", 2);
    parser->addArgument<uint>("rate", 'r', "the frame rate (frames per second)", 1);
    parser->addArgument<double>("step", 's', "the number of steps per frame", 3000);
    parser->addArgument<double>("ba", 'a', "the number of b / a about ecllipse", 1);
    return parser;
}

int main(int argc, char * argv[]) {
    auto parser = BuildArgsParser();
    parser->parse(argc, argv);

    const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
    const auto dim    = std::any_cast<int>(parser->getValueByName("dim"));
    const auto test   = std::any_cast<int>(parser->getValueByName("test"));
    const auto begin  = std::any_cast<uint>(parser->getValueByName("begin"));
    const auto end    = std::any_cast<uint>(parser->getValueByName("end"));
    const auto rate   = std::any_cast<uint>(parser->getValueByName("rate"));
    const auto step   = std::any_cast<double>(parser->getValueByName("step"));
    const auto ba     = std::any_cast<double>(parser->getValueByName("ba"));

    std::unique_ptr<Simulation> smSystem;
    if (dim == 2)
        smSystem = MagneticSandBuilder::build<2>(test, ba);
    else if (dim == 3)
        smSystem = MagneticSandBuilder::build<3>(test, ba);
    else {
        std::cerr << "Error: [main] encountered invalid dimension." << std::endl;
        std::exit(-1);
    }
    auto simulator = std::make_unique<Simulator>(output, begin, end, rate, step, smSystem.get());
    simulator->Simulate();

    return 0;
}
