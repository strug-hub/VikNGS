#include "SimPlotWindow.h"
#include "ui_simplotwindow.h"
#include "TableDisplayWindow.h"

SimPlotWindow::SimPlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SimPlotWindow)
{
    ui->setupUi(this);

    setWindowIcon(QIcon(":icon.svg"));
    connect(ui->simplot_power, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMovePlot1(QMouseEvent*)));
    connect(ui->simplot_power, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickPlot1(QMouseEvent*)));
    connect(ui->simplot_plot2, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(mouseMovePlot2(QMouseEvent*)));
    connect(ui->simplot_plot2, SIGNAL(mouseRelease(QMouseEvent*)), this, SLOT(mouseClickPlot2(QMouseEvent*)));

    QColor blue = QColor(106, 176, 203);
    QColor orange = QColor(255, 187, 63);
    QColor pink = QColor(204, 121, 167);
    QColor green = QColor(106, 150, 71);

    colours.push_back(blue);
    colours.push_back(pink);
    colours.push_back(orange);
    colours.push_back(green);
}

SimPlotWindow::~SimPlotWindow(){
    delete ui;
}

QString timeToString(double time){

    if(time < 120)
        return QString::number(time, 'g', 2) + " sec";

    int t = std::round(time);
    int hour = t/3600; t = t%3600;
    int min = t/60; t = t%60;
    int sec = t;

    QString str = "";
    if(hour > 0)
        str = str + QString::number(hour) + " hr ";
    if(hour > 0 || min > 0)
        str = str + QString::number(min) + " min ";

    str = str + QString::number(sec) + " sec";

    return str;
}

void SimPlotWindow::initialize(Data& results, SimulationRequest& req, QString title){

    this->nsteps = req.steps;
    this->testsPerStep = std::round(1.0 * results.tests.size() / (1.0 * nsteps));
    this->alpha = 0.05;

    if(req.underNull()){
        this->yAxisLabel = "Type I Error";
        ui->simplot_measureGrp->setTitle("Type I Error");
    }
    else{
        this->yAxisLabel = "Power";
        ui->simplot_measureGrp->setTitle("Power");
        ui->simplot_plotTab->setTabEnabled(1, false);
    }

    QString timing = "Simulation time: " + timeToString(results.processingTime) + "\n"
            "Testing time: " + timeToString(results.evaluationTime);

    ui->simplot_timeLbl->setText(timing);

    this->xAxisLabel = "Sample Size";
    this->setWindowTitle("Plotter - " + title);

    getPvalues(results.variants);

    this->result = results;
    this->request = req;

    this->powerIndex = 0;
    this->qqIndex = -1;

    buildPowerPlot();
    buildLegend();

    buildQQPlot(powerIndex);
    updatePowerValues(0);
}

void SimPlotWindow::getPvalues(std::vector<VariantSet>& variants){
    if(variants.size() < 1)
        return;

    for(size_t i = 0; i < variants[0].nPvals(); i++){
        std::vector<double> pval_i(variants.size());

        for(size_t j = 0; j < variants.size(); j++){
            if(variants[j].nPvals() > i)
                pval_i[j] = variants[j].getPval(i);
        }
        std::sort(pval_i.begin(), pval_i.end());
        pvalues.push_back(pval_i);
    }
}

void SimPlotWindow::updateSampleTable(int index){

    if(index < 0) index = 0;

    QTableWidget* table = ui->simplot_sampleTbl;

    for(int i = 0; i < request.groups.size(); i++){

        if(table->rowCount() < i+1)
            table->insertRow(i);

        //table->setItem(i+1, 0, new QTableWidgetItem(QString::number(request.groups[i].index)));
        table->setItem(i+1, 0, new QTableWidgetItem(QString::number(request.groups[i].getSampleSize(index))));
        table->setItem(i+1, 1, new QTableWidgetItem(QString::fromStdString(request.groups[i].getCohort())));

        QString depth = "N(" + QString::number(request.groups[i].meanDepth) + ", " +
                QString::number(request.groups[i].sdDepth) + ")";

        table->setItem(i+1, 2, new QTableWidgetItem(depth));
        table->setItem(i+1, 3, new QTableWidgetItem(QString::number(request.groups[i].errorRate)));

    }
}

void SimPlotWindow::mouseClickPlot1(QMouseEvent *event){
    int closestIndex = findClosestPoint(ui->simplot_power, event);
    if(closestIndex >= 0){
        powerIndex = closestIndex;
        buildQQPlot(closestIndex);
        updatePowerValues(closestIndex);
    }
}

void SimPlotWindow::mouseClickPlot2(QMouseEvent *event){
    if(powerIndex >= 0){
        int closestIndex = findClosestPoint(ui->simplot_plot2, event, true);
        qqIndex = closestIndex;
        buildQQPlot(powerIndex, closestIndex);
    }
}

void SimPlotWindow::mouseMovePlot1(QMouseEvent *event){

    int closestIndex = findClosestPoint(ui->simplot_power, event);
    updateSampleTable(closestIndex);

    if(closestIndex >= 0){

        //remove old highlight layer
        ui->simplot_power->removeGraph(ui->simplot_power->graphCount()-1);
        ui->simplot_power->addGraph();

        QVector<double> x;
        QVector<double> y;

        for(int j = 0; j < ui->simplot_power->graphCount()-1; j++){
            x.push_back(ui->simplot_power->graph(j)->data().data()->at(closestIndex)->key);
            y.push_back(ui->simplot_power->graph(j)->data().data()->at(closestIndex)->value);
        }

        ui->simplot_power->graph()->setData(x, y);
        ui->simplot_power->graph()->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssDisc,
                                    highlight, Qt::white, 12));
        ui->simplot_power->graph()->setLineStyle(QCPGraph::LineStyle::lsNone);

        ui->simplot_power->replot();

        //if(!QToolTip::isVisible())
        //    QToolTip::showText(QPoint(event->globalX(), event->globalY()), QString("hello"));
    }
    else{
        //QToolTip::hideText();
        //remove old highlight layer
        ui->simplot_power->removeGraph(ui->simplot_power->graphCount()-1);
        ui->simplot_power->addGraph();
        ui->simplot_power->replot();
    }
}

void SimPlotWindow::mouseMovePlot2(QMouseEvent *event){

    int closestGraphIndex = findClosestPoint(ui->simplot_plot2, event, true);

    if(qqIndex != closestGraphIndex)
        buildQQPlot(powerIndex, closestGraphIndex);
}

int SimPlotWindow::findClosestPoint(QCustomPlot *plot, QMouseEvent *event, bool getGraphIndex){

    double maxDist = 200;

    double x = event->pos().x();
    double y = event->pos().y();

    double minDist = 9999999999;
    int minIndex;

    for(int j = 0; j < plot->graphCount()-1; j++){
        for(int i = 0; i <  plot->graph(j)->dataCount(); i++){

            double dataCoordx = plot->graph(j)->data().data()->at(i)->key;
            double dataCoordy = plot->graph(j)->data().data()->at(i)->value;

            double dataX = plot->xAxis->coordToPixel(dataCoordx);
            double dataY = plot->yAxis->coordToPixel(dataCoordy);

            double dist = pow(dataX - x, 2) + pow(dataY - y, 2);

            if(dist < minDist){
                minDist = dist;
                if(getGraphIndex)
                    minIndex = j;
                else
                    minIndex = i;
            }
        }
    }
    if(minDist > maxDist)
        return -1;

    return minIndex;
}

double SimPlotWindow::calculatePower(int index, double alpha){
    double success = 0;
    double npval = pvalues[index].size();
    for(size_t i = 0; i < npval; i++)
        if(pvalues[index][i] <= alpha)
            success++;

    return success/npval;
}

void SimPlotWindow::updatePowerValues(int index){

    ui->simplot_powerTbl->setColumnCount(2);
    ui->simplot_powerTbl->setRowCount(testsPerStep);

    for(int i = 0; i < testsPerStep; i++){
        double power = ui->simplot_power->graph(i)->data().data()->at(index)->value;
        ui->simplot_powerTbl->setItem(i, 1, new QTableWidgetItem(QString::number(power)));
        QTableWidgetItem * c = new QTableWidgetItem(QString::fromStdString(result.tests[i].toShortString()));
        c->setBackgroundColor(colours[i]);
        ui->simplot_powerTbl->setItem(i, 0, c);
    }
}

void SimPlotWindow::on_simplot_alphaDial_valueChanged(int value){
    alpha = std::pow(10, -value/1000.0);
    ui->simplot_alphaTxt->setText(QString::number(alpha));
    buildPowerPlot();
    updateAlphaLine();
    updatePowerValues(powerIndex);
}

void SimPlotWindow::on_simplot_alphaTxt_textChanged(const QString &arg1){

    bool ok = false;
    double x = arg1.toDouble(&ok);
    if(ok && x >= 0 && x <=1){
        alpha = x;

        buildPowerPlot();
        updateAlphaLine();
        updatePowerValues(powerIndex);
    }
}

void SimPlotWindow::on_pushButton_clicked(bool checked)
{
    TableDisplayWindow *table = new TableDisplayWindow();
    QString title = "Table";

    std::vector<int> testsToShow;

    for(int i = 0; i < this->testsPerStep; i++)
        testsToShow.push_back(this->powerIndex * this->testsPerStep + i);

    table->initialize(title, &result, testsToShow);
    table->show();
}

void SimPlotWindow::on_simplot_pdfBtn_pressed()
{

    QFileDialog dialog;
    dialog.setFileMode(QFileDialog::Directory);
    dialog.setOption(QFileDialog::ShowDirsOnly);
    QString fileName = QFileDialog::getSaveFileName(this, tr("Export Directory"),
                                                     lastSaveDir, tr("PDF (*.pdf)"));

    if(fileName.size() > 0){
        lastSaveDir = fileName;

        if(ui->simplot_plotTab->currentIndex() == 0)
            saveAsPdf(ui->simplot_power, fileName);
        else if (ui->simplot_plotTab->currentIndex() == 1)
            saveAsPdf(ui->simplot_plot2, fileName, powerIndex);

    }

}

void SimPlotWindow::saveAsPdf(QCustomPlot* plot, QString fileName, int stepIndex){

    //set up PDF
    QPrinter printer;
    printer.setFullPage(true);
    printer.setPaperSize(QPrinter::A4);
    printer.setOrientation(QPrinter::Portrait);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(fileName);

    //print out plot
    double width = printer.width()*0.8;
    double height = (width/plot->width())*plot->height();
    QTextDocument doc;
    QCPDocumentObject *plotObjectHandler = new QCPDocumentObject(this);
    doc.documentLayout()->registerHandler(QCPDocumentObject::PlotTextFormat, plotObjectHandler);

    QTextCursor cursor(&doc);
    cursor.movePosition(QTextCursor::End);
    cursor.insertText(QString(QChar::ObjectReplacementCharacter), QCPDocumentObject::generatePlotFormat(plot, width, height));

    //print out legend
    QVector<QLabel*> labels;
    labels.push_back(ui->simplot_legend1Lbl);
    labels.push_back(ui->simplot_legend2Lbl);
    labels.push_back(ui->simplot_legend3Lbl);
    labels.push_back(ui->simplot_legend4Lbl);
    QVector<QWidget*> colours;
    colours.push_back(ui->simplot_legend1Sqr);
    colours.push_back(ui->simplot_legend2Sqr);
    colours.push_back(ui->simplot_legend3Sqr);
    colours.push_back(ui->simplot_legend4Sqr);

    int legendSize = 0;
    for(int i = 0; i < labels.size(); i++)
        if(!labels[i]->isHidden())
            legendSize++;

    QTextTableFormat tableFormat;
    //tableFormat.setHeaderRowCount(1);
    //tableFormat.setAlignment(Qt::AlignHCenter);
    tableFormat.setCellPadding(10);
    tableFormat.setCellSpacing(0);
    tableFormat.setBorder(1);
    tableFormat.setBorderBrush(QBrush(Qt::SolidPattern));
    tableFormat.clearColumnWidthConstraints();
    QTextTable *textTable = cursor.insertTable(1, 2*legendSize, tableFormat);

    int index = 0;
    for(int i = 0; i < labels.size(); i++){
        if(!labels[i]->isHidden()){
            QTextTableCell cell = textTable->cellAt(0, (i*2)+1);
            QTextCursor cellCursor = cell.firstCursorPosition();
            cellCursor.insertText(labels[i]->text());
            QTextCharFormat cellFormat;
            cellFormat.setBackground(colours[i]->palette().background().color());
            QTextTableCell colourCell = textTable->cellAt(0, i*2);
            colourCell.setFormat(cellFormat);
            index++;
        }
    }

    cursor.movePosition(QTextCursor::End);
    cursor.insertText("\n\n\n");
    //print out table
    QTextTableFormat sampleTableFormat;
    sampleTableFormat.setHeaderRowCount(1);
    //tableFormat.setAlignment(Qt::AlignHCenter);
    sampleTableFormat.setCellPadding(10);
    sampleTableFormat.setCellSpacing(0);
    sampleTableFormat.setBorder(1);
    sampleTableFormat.setBorderBrush(QBrush(Qt::SolidPattern));
    sampleTableFormat.clearColumnWidthConstraints();
    QTextTable *sampleTable = cursor.insertTable(request.groups.size()+1, 6, sampleTableFormat);

    sampleTable->cellAt(0, 0).firstCursorPosition().insertText("Group");
    sampleTable->cellAt(0, 1).firstCursorPosition().insertText("Sample Size");
    sampleTable->cellAt(0, 2).firstCursorPosition().insertText("Cohort");
    sampleTable->cellAt(0, 3).firstCursorPosition().insertText("Mean Read Depth");
    sampleTable->cellAt(0, 4).firstCursorPosition().insertText("Read Depth SD");
    sampleTable->cellAt(0, 5).firstCursorPosition().insertText("Error Rate");

    QTextCharFormat sampleTableHeaderFormat;
    sampleTableHeaderFormat.setBackground(QColor("#DADADA"));
    for (size_t i = 0; i < 6; i++)
        sampleTable->cellAt(0, i).setFormat(sampleTableHeaderFormat);

    for(size_t i = 0; i < request.groups.size(); i++){
        sampleTable->cellAt(i+1, 0).firstCursorPosition().insertText(
                    QString::number(request.groups[i].index));
        QString sampleSize = QString::number(request.groups[i].getSampleSize(0));
        if(request.groups[i].getSampleSize(0) != request.groups[i].getSampleSize(request.steps))
            sampleSize = sampleSize + ":" + QString::number(request.groups[i].getSampleSize(request.steps));

        if(stepIndex >= 0)
            sampleSize = QString::number(request.groups[i].getSampleSize(stepIndex));

        sampleTable->cellAt(i+1, 1).firstCursorPosition().insertText(
                    (sampleSize));
        sampleTable->cellAt(i+1, 2).firstCursorPosition().insertText(
                    QString::fromStdString(request.groups[i].getCohort()));
        sampleTable->cellAt(i+1, 3).firstCursorPosition().insertText(
                    QString::number(request.groups[i].meanDepth));
        sampleTable->cellAt(i+1, 4).firstCursorPosition().insertText(
                    QString::number(request.groups[i].sdDepth));
        sampleTable->cellAt(i+1, 5).firstCursorPosition().insertText(
                    QString::number(request.groups[i].errorRate));
    }


    doc.print(&printer);
}
