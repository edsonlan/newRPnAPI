#ifndef _OBSERVER_
#define _OBSERVER_

class Subject;

class Observer {
    private:
    protected:
    public:
        virtual ~Observer();

        virtual void change() = 0;
};

#endif // _OBSERVER_

