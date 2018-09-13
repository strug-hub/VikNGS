#include "PlotWindow.h"
#include "ui_plotwindow.h"

PlotWindow::PlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PlotWindow)
{
    ui->setupUi(this);
    qApp->installEventFilter(this);

    setWindowIcon(QIcon(":icon.svg"));

    focusedVar = -1;

    connect(ui->plot_chrPlt, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMoveChromosome(QMouseEvent*)));
    connect(ui->plot_genomePlt, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMoveGenome(QMouseEvent*)));
    connect(ui->plot_genomePlt, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickGenome(QMouseEvent*)));
    connect(ui->plot_chrPlt, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickChromosome(QMouseEvent*)));

    connect(ui->plot_chrPlt->xAxis, SIGNAL(rangeChanged(const QCPRange&)), this, SLOT(updateChromosomeRange(const QCPRange&)));

    ui->plot_chrPlt->plotLayout()->insertRow(0);

}

PlotWindow::~PlotWindow()
{
    delete ui;
}

void PlotWindow::initialize(Data& result, QString title){

    this->result = result;

    ui->plot_title->setText(title);
    this->setWindowTitle("Plotter - " + title);
    createChromosomes(this->result.variants);
    focusedChr = chrNames[0];
    buildGenomePlot();
    buildChromosomePlot(focusedChr);

    std::vector<int> testsToShow;
    for(size_t i = 0; i < this->result.tests.size(); i++) testsToShow.push_back(i);

    QString table_title = "Table";
    tableView->initialize(table_title, &this->result, testsToShow);
}

void PlotWindow::createChromosomes(std::vector<VariantSet>& variantSets){

    for(size_t i = 0; i < variantSets.size(); i++){
        if(variantSets[i].nPvals() < 1) continue;
        for(size_t j = 0; j < variantSets[i].size(); j++){
            Variant* v = &variantSets[i].getVariants()->at(j);
            if(!v->isValid()) continue;
            QString chr = QString::fromStdString(v->getChromosome());
            if(chromosomes.contains(chr))
                chromosomes[chr].addVariant(v, &variantSets[i]);
            else{
                Chromosome newChr(chr, v, &variantSets[i], grey1);
                chromosomes[chr] = newChr;
                chrNames.push_back(chr);
            }
        }
    }

    std::sort(chrNames.begin(), chrNames.end(), chrCompare);
}


void PlotWindow::updateVariantInfo(int variantIndex){

   if(variantIndex < 0)
       return;

   Variant* variant = chromosomes[focusedChr].getVariant(variantIndex);

   QString pval = QString::number(std::pow(10, -chromosomes[focusedChr].getPval(variantIndex)));
   QString chr = QString::fromStdString(variant->getChromosome());
   QString pos = QString::number(variant->getPosition());
   QString ref = QString::fromStdString(variant->getRef());
   QString alt = QString::fromStdString(variant->getAlt());

   if(ui->plot_refTxt->text() != ref)
        ui->plot_refTxt->setText(ref);
   if(ui->plot_altTxt->text() != alt)
        ui->plot_altTxt->setText(alt);

    QString info = chr + " - " + pos + ": " + pval;
    if(ui->plot_variantInfoLbl->text() != info)
        ui->plot_variantInfoLbl->setText(info);
}

void PlotWindow::moveRectangle(QCPItemRect* rect, QString chr, double lower, double upper){

    double min = chromosomes[chr].getMinPos();
    double max = chromosomes[chr].getMaxPos();

    double offset = chromosomes[chr].getOffset();
    double top = ui->plot_genomePlt->yAxis->range().upper;

    if(lower >= 0)
        rect->topLeft->setCoords(offset + (lower - min), top);
    else
        rect->topLeft->setCoords(offset, top);

    if(upper >= 0)
        rect->bottomRight->setCoords(offset + (max - min) - (max - upper), 0);
    else
        rect->bottomRight->setCoords(offset + (max - min), 0);

}

QCPGraph* PlotWindow::getGraphByName(QCustomPlot *plot, QString name){
    for (int i=0; i < plot->graphCount(); ++i)
      if (plot->graph(i)->name() == name)
        return plot->graph(i);

    return plot->graph();
}

void PlotWindow::resetColor(QString chrName){
    getGraphByName(ui->plot_genomePlt, chrName)->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                chromosomes[chrName].getColour(), Qt::white, 2));
}

QString PlotWindow::getChromUnderCursor(QMouseEvent *event){

    double x = ui->plot_genomePlt->xAxis->pixelToCoord(event->pos().x());
    double y = ui->plot_genomePlt->yAxis->pixelToCoord(event->pos().y());


    if(y > 0 && x > 0 &&
            y < ui->plot_genomePlt->yAxis->range().upper &&
            x < ui->plot_genomePlt->xAxis->range().upper){
        double offset = 0;

        for(int i = 0; i < chrNames.size(); i++){
            offset += chromosomes[chrNames[i]].getSpan();
            if (offset > x)
                return chromosomes[chrNames[i]].getName();
            }
        }
    return "";
}


void PlotWindow::on_plot_genotypeBtn_pressed()
{
    if(!tableView->isVisible())
        tableView->show();
}
