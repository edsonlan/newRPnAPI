#ifndef _SUBJECT_
#define _SUBJECT_

#include <set>
#include "Observer.h"

class Subject {
    private:
    protected:
        std::set<Observer*> observer_;
    public:
        virtual ~Subject();

        virtual void add(Observer *obs){
            observer_.insert(obs);
            return;
        }

        virtual void remove(Observer *obs){
            observer_.erase(observer_.find(obs));
            return;
        }

        virtual void notify(){
             std::set<Observer*>::iterator it;
            for (it = observer_.begin(); it != observer_.end(); it++) (*it)->change();

            return;
        }
};

#endif // _SUBJECT_

