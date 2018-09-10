#include "PlotWindow.h"
#include "ui_plotwindow.h"

PlotWindow::PlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PlotWindow)
{
    ui->setupUi(this);
    qApp->installEventFilter(this);

    setWindowIcon(QIcon(":icon.svg"));

    VariantSet null;
    nullVariant = null;
    focusedVar = &nullVariant;

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

    ui->plot_title->setText(title);
    this->setWindowTitle("Plotter - " + title);
    createChromosomes(result.variants);
    focusedChr = chrNames[0];
    buildGenomePlot();
    buildChromosomePlot(focusedChr);

    this->result = result;

    QString table_title = "Table";
    tableView->initialize(table_title, &this->result);
}

void PlotWindow::createChromosomes(std::vector<VariantSet>& variantSets){

    for(size_t i = 0; i < variantSets.size(); i++){
        VariantSet* vs = &variantSets[i];
        QString chr = QString::fromStdString(vs->getChromosome());
        if(chromosomes.contains(chr))
            chromosomes[chr].addVariant(vs);
        else{
            Chromosome newChr(chr, vs, grey1);
            chromosomes[chr] = newChr;
            chrNames.push_back(chr);
        }
    }

    std::sort(chrNames.begin(), chrNames.end(), chrCompare);
}


void PlotWindow::updateVariantInfo(VariantSet* variant){

   QString chr = QString::fromStdString(variant->getChromosome());
   QString pos = QString::number(variant->getMinPos());
   QString ref = QString::fromStdString("TODO");
   QString alt = QString::fromStdString("TODO");


   if(ui->plot_refTxt->text() != ref)
        ui->plot_refTxt->setText(ref);
   if(ui->plot_altTxt->text() != alt)
        ui->plot_altTxt->setText(alt);

    QString info = chr + " - " + pos + ":";
    if(ui->plot_variantInfoLbl->text() != info)
        ui->plot_variantInfoLbl->setText(info);
}

void PlotWindow::moveRectangle(QCPItemRect *rect, QString chrName, double lower, double upper){

    double width = chromosomes[chrName].getSpan();
    double start = chromosomes[chrName].getOffset();
    double top = ui->plot_genomePlt->yAxis->range().upper;

    double lowerShift = 0;
    if(lower > 0)
        lowerShift = lower;

    rect->topLeft->setCoords( start + lowerShift, top );

    double upperShift = 0;
    if(upper > 0)
        upperShift = std::max(0.0, chromosomes[chrName].getMaxPos() - upper);

    rect->bottomRight->setCoords(start + width - upperShift, 0);
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
