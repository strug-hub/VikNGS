#ifndef CHECKTREE_H
#define CHECKTREE_H

#include <QTreeWidgetItem>

class CheckTree : public QObject
{
    Q_OBJECT

public:
    CheckTree();
    ~CheckTree(){}

private:
    std::vector<CheckTree*> children;
    QString name;
    int index;
    bool checked = true;
    QTreeWidgetItem* treeItem;

public:

    void initalize(QString name, int index=-1, bool checked=true){
        //this->name=name;
        this->index=index;
        this->checked=checked;
    }
    void addChild(QString childName, int childIndex=-1, bool childChecked = true){
        CheckTree* newChild;
        newChild->initalize(childName, childIndex, childChecked);
        children.push_back(newChild);
    }

    QTreeWidgetItem* draw(){

        treeItem = new QTreeWidgetItem();

        for(CheckTree* child : children)
           treeItem->addChild(child->draw());

        treeItem->setCheckState(0,Qt::Checked);
        treeItem->setText(0, name);

        treeItem->setExpanded(true);
        return treeItem;
    }


public slots:
    bool updateState(){
        for(CheckTree* child : children)
            child->updateState();

        bool allChildrenChecked=true;
        bool someChildrenChecked=false;

        for(CheckTree* child : children){
            allChildrenChecked=allChildrenChecked && child->checked;
            someChildrenChecked=someChildrenChecked || child->checked;
        }

        if(someChildrenChecked){
            checked=true;

           if(allChildrenChecked)
               treeItem->setCheckState(0, Qt::Checked);
           else
               treeItem->setCheckState(0, Qt::PartiallyChecked);
        }
        else{
            checked=false;
            treeItem->setCheckState(0, Qt::Unchecked);
        }

        return false;
    }
} ;


#endif // CHECKTREE_H
