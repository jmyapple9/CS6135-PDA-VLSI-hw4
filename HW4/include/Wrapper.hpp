#pragma once
#include <memory>

// Wrapped by Jinkela 2023.11.22

class Module;
class Net;
class Pin;
class Placement;
class Rectangle;
class Row;

namespace wrapper
{
    class Module
    {
    public:
        Module(::Module &_origin);

        /////////////////////////////////////////////
        // get
        /////////////////////////////////////////////
        const char *name();
        double x();
        double y();
        double width();
        double height();
        bool isFixed();

        double centerX();
        double centerY();
        double area();

        Rectangle rectangle();

        const char *orientString();

        /////////////////////////////////////////////
        // set
        /////////////////////////////////////////////
        void setPosition(double x, double y);
        void setCenterPosition(double x, double y);

        /////////////////////////////////////////////
        // get (for pins of this modules)
        /////////////////////////////////////////////
        unsigned numPins();
        Pin &pin(unsigned index);

        ::Module &getOrigin();
        const ::Module &getOrigin() const;

    private:
        ::Module &_origin;
    };

    class Placement
    {
    public:
        Placement();
        ~Placement();
        /////////////////////////////////////////////
        // input/output
        /////////////////////////////////////////////
        void readBookshelfFormat(const char *const filePathName, const char *const plFileName);
        void outputBookshelfFormat(const char *const filePathName); // output file function

        /////////////////////////////////////////////
        // get
        /////////////////////////////////////////////
        const char *name();
        const char *plname();
        double computeHpwl();
        double computeTotalNetLength(int cellid);
        Rectangle rectangleChip();
        double boundryTop();
        double boundryLeft();
        double boundryBottom();
        double boundryRight();

        /////////////////////////////////////////////
        // operation
        /////////////////////////////////////////////
        void moveDesignCenter(double xOffset, double yOffset);

        /////////////////////////////////////////////
        // get design objects/properties
        /////////////////////////////////////////////
        Module module(unsigned moduleId);
        Net &net(unsigned netId);
        Pin &pin(unsigned pinId);
        Row &row(unsigned rowId);

        double getRowHeight();

        unsigned numModules();
        unsigned numNets();
        unsigned numPins();
        unsigned numRows();

        ::Placement &getOrigin();
        const ::Placement &getOrigin() const;

    private:
        std::unique_ptr<::Placement> _origin;
    };
}