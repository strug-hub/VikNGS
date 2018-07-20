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

    connect(ui->plot_chrPlt->xAxis, SIGNAL(rangeChanged(const QCPRange&)), this, SLOT(updateChromosomeRange(const QCPRange&)));

    ui->plot_chrPlt->plotLayout()->insertRow(0);

}

PlotWindow::~PlotWindow()
{
    delete ui;
}


void PlotWindow::initialize(Result& result, QString title){

    ui->plot_title->setText(title);
    this->setWindowTitle("Plotter - " + title);
    createChromosomes(result.variants);
    //setRandomChromosomes();
    focusedChr = chrNames[0];
    buildGenomePlot();
    buildChromosomePlot(focusedChr);

    this->result = result;

    QString table_title = "Table";
    tableView->initialize(table_title, &this->result);
}

void PlotWindow::initialize(int nvariants, QString title){

    ui->plot_title->setText(title);
    this->setWindowTitle("Plotter - " + title);
    setRandomChromosomes(nvariants);
    focusedChr = chrNames[0];
    buildGenomePlot();
    buildChromosomePlot(focusedChr);

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
        p.push_back(-log10(variant.getPval(0)));
    }
    if(focusedVar.isValid()){
        l.push_back(focusedVar.pos);
        p.push_back(-log10(focusedVar.getPval(0)));
    }

    ui->plot_chrPlt->graph()->setData(l,p);
    ui->plot_chrPlt->replot();
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

void PlotWindow::createChromosomes(std::vector<Variant>& variants){

    for(Variant v : variants){
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
        chromosomes[chrNames[i]].setOffset(offset);
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

    focusRect = new QCPItemRect(ui->plot_genomePlt);
    focusRect->setLayer("overlay");
    focusRect->topLeft->setType(QCPItemPosition::ptPlotCoords);
    focusRect->topLeft->setAxisRect( ui->plot_genomePlt->axisRect() );
    focusRect->bottomRight->setType(QCPItemPosition::ptPlotCoords);
    focusRect->bottomRight->setAxisRect( ui->plot_genomePlt->axisRect() );
    focusRect->setPen(QPen(focus));
    QColor transFocus = focus;
    transFocus.setAlpha(50);
    focusRect->setBrush(QBrush(transFocus));

    zoomRect = new QCPItemRect(ui->plot_genomePlt);
    zoomRect->setLayer("overlay");
    zoomRect->topLeft->setType(QCPItemPosition::ptPlotCoords);
    zoomRect->topLeft->setAxisRect( ui->plot_genomePlt->axisRect() );
    zoomRect->bottomRight->setType(QCPItemPosition::ptPlotCoords);
    zoomRect->bottomRight->setAxisRect( ui->plot_genomePlt->axisRect() );
    zoomRect->setPen(Qt::NoPen);
    zoomRect->setBrush(QBrush(transFocus));

    ui->plot_genomePlt->yAxis->setLabel("-log(p)");

    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(tickValues, tickLabels);
    ui->plot_genomePlt->xAxis->setTicker(textTicker);
    ui->plot_genomePlt->xAxis->setTickLabelRotation(90);
    ui->plot_genomePlt->xAxis->grid()->setPen(Qt::NoPen);

    ui->plot_genomePlt->rescaleAxes();
    ui->plot_genomePlt->yAxis->setRangeLower(0);

    double ymax = ui->plot_genomePlt->yAxis->range().upper;
    if(ymax < 7.5){
        ui->plot_genomePlt->yAxis->setRangeUpper(7.5);
        ymax = 7.5;
    }
    ui->plot_genomePlt->yAxis->setRangeUpper(ymax + ymax * 0.05);

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

    moveRectangle(focusRect, ui->plot_genomePlt->graph(0)->name());
    moveRectangle(zoomRect, ui->plot_genomePlt->graph(0)->name());

    ui->plot_genomePlt->replot();
}


void PlotWindow::buildChromosomePlot(QString chrName){
    ui->plot_chrPlt->clearPlottables();

    focusedVar = nullVariant;
    focusedChr = chrName;
    Chromosome chr = chromosomes[chrName];

    ui->plot_chrPlt->addGraph();

    QVector<double> pos = chr.getPositions();
    QVector<double> pvals = chr.getPvals();

    QVector<double> pos_range;
    QVector<double> pvals_range;

    for(int i =0; i<pos.size(); i++){
        pos_range.push_back(pos[i]);
        pvals_range.push_back(pvals[i]);
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
    ui->plot_chrPlt->yAxis->setRangeLower(0);

    double ymax = ui->plot_chrPlt->yAxis->range().upper;
    if(ymax < 7.5){
        ui->plot_chrPlt->yAxis->setRangeUpper(7.5);
        ymax = 7.5;
    }
    ui->plot_chrPlt->yAxis->setRangeUpper(ymax + ymax * 0.05);

    QPen hpen = QPen(Qt::DashDotLine);
    hpen.setColor( QColor::fromRgb(210, 80, 80));
    horizontal = new QCPItemLine(ui->plot_chrPlt);
    horizontal->setPen(hpen);
    horizontal->start->setCoords(0,1.301);
    horizontal->end->setCoords(chr.getMaxPos()*2,1.301);

    QPen bpen = QPen(Qt::DashDotLine);
    bpen.setColor( QColor::fromRgb(250, 145, 145));
    QCPItemLine *bonferroni = new QCPItemLine(ui->plot_chrPlt);
    bonferroni->setPen(bpen);
    bonferroni->start->setCoords(0,7.301);
    bonferroni->end->setCoords(chr.getMaxPos()*2, 7.301);

    ui->plot_chrPlt->addGraph();
    ui->plot_chrPlt->graph()->setName("annotationLayer");
    ui->plot_chrPlt->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);
    ui->plot_chrPlt->graph()->setScatterStyle(
                QCPScatterStyle(QCPScatterStyle::ssDisc,
                                highlight, Qt::white, 6));

    updateVariantHighlightLayer(nullVariant);

    ui->plot_chrPlt->setInteraction(QCP::iRangeDrag, true);
    ui->plot_chrPlt->setInteraction(QCP::iRangeZoom, true);
    ui->plot_chrPlt->axisRect()->setRangeDrag(Qt::Horizontal);
    ui->plot_chrPlt->axisRect()->setRangeZoom(Qt::Horizontal);


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
        s.addPval((float) rand()/RAND_MAX, Test::NONE);


        c.addVariant(s);
    }

    return c;
}

void PlotWindow::setRandomChromosomes(int n){

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

void PlotWindow::on_plot_genotypeBtn_pressed()
{
    if(!tableView->isVisible())
        tableView->show();
}

