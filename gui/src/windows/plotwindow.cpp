#include "plotwindow.h"
#include "ui_plotwindow.h"

PlotWindow::PlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PlotWindow)
{
    ui->setupUi(this);
    qApp->installEventFilter(this);

    setWindowIcon(QIcon(":icon.svg"));

    nullVariant.setInvalid("null");
    focusedVar = nullVariant;
    //connect(ui->plot_window, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMoveWindow(QMouseEvent*)));

    connect(ui->plot_genomePlt, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMoveGenome(QMouseEvent*)));
    connect(ui->plot_genomePlt, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickGenome(QMouseEvent*)));
    connect(ui->plot_chrPlt, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickChromosome(QMouseEvent*)));
    connect(ui->plot_chrPlt, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMoveChromosome(QMouseEvent*)));
    connect(ui->plot_chrPlt, SIGNAL(mouseWheel(QWheelEvent*)), this, SLOT(mouseScrollChromosome(QWheelEvent*)));
    connect(ui->plot_zoomBar, SIGNAL(rangeChanged(double,double)), this, SLOT(zoomChromosomePlot(double,double)));
    ui->plot_chrPlt->plotLayout()->insertRow(0);
}

PlotWindow::~PlotWindow()
{
    delete ui;
}

void PlotWindow::zoomChromosomePlot(double min, double max){
    buildChromosomePlot(focusedChr, min, max);
}

void PlotWindow::initialize(QVector<Variant> variants, QString title){

    ui->plot_title->setText(title);
    this->setWindowTitle("Plotter - " + title);
    createChromosomes(variants);
    //setRandomChromosomes();
    focusedChr = chrNames[0];
    buildGenomePlot();
    buildChromosomePlot(focusedChr, -1, -1);

}

void PlotWindow::updateVariantInfo(Variant variant){

   if(!variant.isValid()){
       ui->plot_variantInfoLbl->setText("");
       ui->plot_altTxt->setText("");
       ui->plot_refTxt->setText("");

       return;
    }

   QString chr = QString::fromStdString(variant.chr);
   QString pos = QString::number(variant.pos);
   QString ref = QString::fromStdString(variant.ref);
   QString alt = QString::fromStdString(variant.alt);


   if(ui->plot_refTxt->text() != ref)
        ui->plot_refTxt->setText(ref);
   if(ui->plot_altTxt->text() != alt)
        ui->plot_altTxt->setText(alt);

    QString info = chr + " - " + pos + ":";
    if(ui->plot_variantInfoLbl->text() != info)
        ui->plot_variantInfoLbl->setText(info);
}

void PlotWindow::updateVariantHighlightLayer(Variant variant){
    QVector<double> l;
    QVector<double> p;

    if(variant.isValid()){
        l.push_back(variant.pos);
        p.push_back(-log10(variant.pvalue));
    }
    if(focusedVar.isValid()){
        l.push_back(focusedVar.pos);
        p.push_back(-log10(focusedVar.pvalue));
    }

    ui->plot_chrPlt->graph()->setData(l,p);
    ui->plot_chrPlt->replot();
}



QCPGraph* PlotWindow::getGraphByName(QCustomPlot *plot, QString name){
    for (int i=0; i < plot->graphCount(); ++i)
      if (plot->graph(i)->name() == name)
        return plot->graph(i);

    return plot->graph();
}

QString PlotWindow::getChromUnderCursor(QMouseEvent *event){

    double x = ui->plot_genomePlt->xAxis->pixelToCoord(event->pos().x());
    double y = ui->plot_genomePlt->yAxis->pixelToCoord(event->pos().y());


    if(y > 0 && x > 0 && y < ui->plot_genomePlt->yAxis->atTop){
        double offset = 0;

        for(int i = 0; i < chrNames.size(); i++){
            offset += chromosomes[chrNames[i]].getSpan();
            if (offset > x)
                return chromosomes[chrNames[i]].getName();
            }
        }
    return "";
}

void PlotWindow::createChromosomes(QVector<Variant> variants){

    QVectorIterator<Variant> it(variants);
    while (it.hasNext()){
        Variant v = it.next();
        QString chr = QString::fromStdString(v.chr);
        if(chromosomes.contains(chr))
            chromosomes[chr].addVariant(v);
        else{
            Chromosome newChr(chr, v, grey1);
            chromosomes[chr] = newChr;
            chrNames.push_back(chr);
        }
    }

    std::sort(chrNames.begin(), chrNames.end(), chrCompare);
}


void PlotWindow::buildGenomePlot(){
    ui->plot_genomePlt->clearPlottables();

    QVector<double> tickValues;
    QVector<QString> tickLabels;

    bool colour1 = true;    

    double offset = 0;

    for( int i = 0; i<chrNames.size(); i++ )
    {
        Chromosome chr = chromosomes[chrNames[i]];

        ui->plot_genomePlt->addGraph();
        ui->plot_genomePlt->graph()->setData(chr.getPositions(offset), chr.getPvals());
        ui->plot_genomePlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
        ui->plot_genomePlt->graph()->setName(chr.getName());

        tickLabels.push_back(chr.getName());
        tickValues.push_back(offset + chr.getSpan()/2);

        QColor toUse = grey1;
        if(!colour1)
            toUse=grey2;

        colour1 = !colour1;

        chromosomes[chrNames[i]].setColour(toUse);

        if(chrNames[i] == focusedChr)
            toUse = focus;

        ui->plot_genomePlt->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    toUse, Qt::white, 2));
        offset += chr.getSpan();       
    }

    ui->plot_genomePlt->yAxis->setLabel("-log(p)");

    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(tickValues, tickLabels);
    ui->plot_genomePlt->xAxis->setTicker(textTicker);
    ui->plot_genomePlt->xAxis->setTickLabelRotation(90);
    ui->plot_genomePlt->xAxis->grid()->setPen(Qt::NoPen);

    ui->plot_genomePlt->rescaleAxes();
    if(ui->plot_genomePlt->yAxis->atTop < 7.5)
        ui->plot_genomePlt->yAxis->setRange(0, 7.5);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(ui->plot_genomePlt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(ui->plot_genomePlt->xAxis->range().upper,1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(ui->plot_genomePlt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(ui->plot_genomePlt->xAxis->range().upper,7.301);

    ui->plot_genomePlt->replot();

}

void PlotWindow::buildChromosomePlot(QString chrName, double min, double max){
    ui->plot_chrPlt->clearPlottables();

    focusedChr = chrName;
    Chromosome chr = chromosomes[chrName];

    ui->plot_chrPlt->addGraph();

    QVector<double> pos = chr.getPositions();
    QVector<double> pvals = chr.getPvals();

    QVector<double> pos_range;
    QVector<double> pvals_range;

    for(int i =0; i<pos.size(); i++){

        if(min < 0 || pos[i] > min)
            if(max < 0 || pos[i] < max){
                pos_range.push_back(pos[i]);
                pvals_range.push_back(pvals[i]);
            }
    }

    ui->plot_chrPlt->graph()->setData(pos_range,pvals_range);
    ui->plot_chrPlt->graph()->setName(chrGraphName);
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

    double pointSize = 1.0/(sqrt(pos_range.size())/100.0);
    if(pointSize > 6) pointSize = 6;
    if(pointSize < 2) pointSize = 2;

    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                focus, Qt::white, pointSize));

    ui->plot_chrPlt->yAxis->setLabel("-log(p)");
    ui->plot_chrPlt->xAxis->setLabel("Position");

    QCPTextElement *title = new QCPTextElement(ui->plot_chrPlt, "Chromosome " + chrName, QFont("sans", 17, QFont::Bold));
    ui->plot_chrPlt->plotLayout()->removeAt(0);
    ui->plot_chrPlt->plotLayout()->addElement(0,0,title);

    ui->plot_chrPlt->rescaleAxes();
    if(ui->plot_chrPlt->yAxis->atTop < 7.5)
        ui->plot_chrPlt->yAxis->setRange(0, 7.5);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(ui->plot_chrPlt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(chr.getMaxPos(),1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(ui->plot_chrPlt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(chr.getMaxPos(), 7.301);

    ui->plot_chrPlt->addGraph();
    ui->plot_chrPlt->graph()->setName("annotationLayer");
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                highlight, Qt::white, 6));

    updateVariantHighlightLayer(nullVariant);

    if(min < 0)
        ui->plot_zoomBar->updateMin(ui->plot_chrPlt->xAxis->range().lower);
    else
        ui->plot_zoomBar->updateZoomMin(min);

    if(max < 0)
        ui->plot_zoomBar->updateMax(ui->plot_chrPlt->xAxis->range().upper);
    else
        ui->plot_zoomBar->updateZoomMax(max);


     ui->plot_chrPlt->replot();
}


//----------------------------------------------------------------------------
//for development purposes
//----------------------------------------------------------------------------

#include <time.h>

Chromosome PlotWindow::generateRandomChromosome(int n, std::string chrom, int maxPos){

    Chromosome c;
    c.setName(QString::fromStdString(chrom));
    c.setColour(QColor::fromRgb(round((float) rand()/RAND_MAX*254),
                                round((float) rand()/RAND_MAX*254),
                                round((float) rand()/RAND_MAX*254)));


    for(int sample = 0; sample < n; sample++){

        Variant s;

        double alt = (double) rand()/RAND_MAX;
        if(alt < 0.25)
            s.alt = "T";
        else if(alt < 0.5)
            s.alt = "C";
        else if(alt < 0.75)
            s.alt = "A";
        else
            s.alt = "G";

        double ref = (double) rand()/RAND_MAX;
        if(ref < 0.25)
            s.ref = "T";
        else if(ref < 0.5)
            s.ref = "C";
        else if(ref < 0.75)
            s.ref = "A";
        else
            s.ref = "G";

        s.pos = (rand() % static_cast<int>(maxPos + 1));
        s.chr = chrom;
        s.pvalue = (float) rand()/RAND_MAX;


        c.addVariant(s);
    }

    return c;
}

void PlotWindow::setRandomChromosomes(){

    int n=5000;
    srand( (unsigned)time( NULL ) );

    chromosomes["chr1"] = generateRandomChromosome(n,"chr1",248956422);
    chromosomes["chr2"] = generateRandomChromosome(n,"chr2",242193529);
    chromosomes["chr3"] = generateRandomChromosome(n,"chr3",198295559);
    chromosomes["chr4"] = generateRandomChromosome(n,"chr4",190214555);
    chromosomes["chr5"] = generateRandomChromosome(n,"chr5",181538259);
    chromosomes["chr6"] = generateRandomChromosome(n,"chr6",170805979);
    chromosomes["chr7"] = generateRandomChromosome(n,"chr7",159345973);
    chromosomes["chr8"] = generateRandomChromosome(n,"chr8",145138636);
    chromosomes["chr9"] = generateRandomChromosome(n,"chr9",138394717);
    chromosomes["chr10"] = generateRandomChromosome(n,"chr10",133797422);
    chromosomes["chr11"] = generateRandomChromosome(n,"chr11",135086622);
    chromosomes["chr12"] = generateRandomChromosome(n,"chr12",133275309);
    chromosomes["chr13"] = generateRandomChromosome(n,"chr13",114364328);
    chromosomes["chr14"] = generateRandomChromosome(n,"chr14",107043718);
    chromosomes["chr15"] = generateRandomChromosome(n,"chr15",101991189);
    chromosomes["chr16"] = generateRandomChromosome(n,"chr16",90338345);
    chromosomes["chr17"] = generateRandomChromosome(n,"chr17",83257441);
    chromosomes["chr18"] = generateRandomChromosome(n,"chr18",80373285);
    chromosomes["chr19"] = generateRandomChromosome(n,"chr19",58617616);
    chromosomes["chr20"] = generateRandomChromosome(n,"chr20",64444167);
    chromosomes["chr21"] = generateRandomChromosome(n,"chr21",46709983);
    chromosomes["chr22"] = generateRandomChromosome(n,"chr22",50818468);
    chromosomes["chrX"] = generateRandomChromosome(n,"chrX",156040895);
    chromosomes["chrY"] = generateRandomChromosome(n,"chrY",57227415);

    chrNames.push_back("chr2");
    chrNames.push_back("chr20");
    chrNames.push_back("chr19");
    chrNames.push_back("chr18");
    chrNames.push_back("chr17");
    chrNames.push_back("chrX");
    chrNames.push_back("chrY");
    chrNames.push_back("chr16");
    chrNames.push_back("chr15");
    chrNames.push_back("chr1");
    chrNames.push_back("chr14");
    chrNames.push_back("chr13");
    chrNames.push_back("chr10");
    chrNames.push_back("chr11");
    chrNames.push_back("chr12");
    chrNames.push_back("chr21");
    chrNames.push_back("chr22");
    chrNames.push_back("chr9");
    chrNames.push_back("chr5");
    chrNames.push_back("chr6");
    chrNames.push_back("chr7");
    chrNames.push_back("chr8");
    chrNames.push_back("chr4");
    chrNames.push_back("chr3");

    qSort(chrNames.begin(), chrNames.end(), chrCompare);


}
