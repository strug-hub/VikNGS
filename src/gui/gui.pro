TEMPLATE = subdirs

#-------------------------------------------------
#
# Project created by QtCreator 2017-08-27T15:59:42
#
#-------------------------------------------------

QT += core gui
QT += printsupport
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = GUI
TEMPLATE = app
CONFIG += embed_manifest_exe
CONFIG += optimize_full

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

RESOURCES += resources/resources.qrc

SOURCES += \
    src/main.cpp \
    src/widgets/qzoombar.cpp \
    src/widgets/qcustomplot.cpp \
    ../Parser/MemoryMapped/MemoryMapped.cpp \
    ../Request.cpp \
    src/windows/MainTab.cpp \
    src/windows/SimulationTab.cpp \
    ../vikNGS.cpp \
    ../Parser/SampleParser.cpp \
    ../Parser/Filter.cpp \
    ../Parser/InputProcess.cpp \
    ../Parser/VariantParser.cpp \
    ../Math/RandomHelper.cpp \
    ../Math/StatisticsHelper.cpp \
    ../Math/VectorHelper.cpp \
    ../Parser/StringTools.cpp \
    ../Math/GeneticsHelper.cpp \
    ../Parser/BEDParser.cpp \
    ../Test/Test.cpp \
    src/windows/PlotWindowPlotter.cpp \
    src/simulation/Simulation.cpp \
    src/windows/SimPlotWindowPlotter.cpp \
    src/windows/TableDisplayWindow.cpp \
    src/windows/TableDisplayWindowColumns.cpp \
    src/simulation/SimulationHelper.cpp \
    src/windows/SimPlotWindow.cpp \
    src/windows/PlotWindow.cpp \
    ../Global.cpp \
    src/log/Log.cpp \
    src/windows/PlotWindowMouse.cpp \
    src/windows/MainWindow.cpp \
    src/widgets/qcpdocumentobject.cpp \
    ../Test/ScoreTestFunctions.cpp

HEADERS += \
    ../Eigen/src/Cholesky/LDLT.h \
    ../Eigen/src/Cholesky/LLT.h \
    ../Eigen/src/Cholesky/LLT_LAPACKE.h \
    ../Eigen/src/Core/arch/AltiVec/Complex.h \
    ../Eigen/src/Core/arch/AltiVec/MathFunctions.h \
    ../Eigen/src/Core/arch/AltiVec/PacketMath.h \
    ../Eigen/src/Core/arch/AVX/Complex.h \
    ../Eigen/src/Core/arch/AVX/MathFunctions.h \
    ../Eigen/src/Core/arch/AVX/PacketMath.h \
    ../Eigen/src/Core/arch/AVX/TypeCasting.h \
    ../Eigen/src/Core/arch/AVX512/MathFunctions.h \
    ../Eigen/src/Core/arch/AVX512/PacketMath.h \
    ../Eigen/src/Core/arch/CUDA/Complex.h \
    ../Eigen/src/Core/arch/CUDA/Half.h \
    ../Eigen/src/Core/arch/CUDA/MathFunctions.h \
    ../Eigen/src/Core/arch/CUDA/PacketMath.h \
    ../Eigen/src/Core/arch/CUDA/PacketMathHalf.h \
    ../Eigen/src/Core/arch/CUDA/TypeCasting.h \
    ../Eigen/src/Core/arch/Default/Settings.h \
    ../Eigen/src/Core/arch/NEON/Complex.h \
    ../Eigen/src/Core/arch/NEON/MathFunctions.h \
    ../Eigen/src/Core/arch/NEON/PacketMath.h \
    ../Eigen/src/Core/arch/SSE/Complex.h \
    ../Eigen/src/Core/arch/SSE/MathFunctions.h \
    ../Eigen/src/Core/arch/SSE/PacketMath.h \
    ../Eigen/src/Core/arch/SSE/TypeCasting.h \
    ../Eigen/src/Core/arch/ZVector/Complex.h \
    ../Eigen/src/Core/arch/ZVector/MathFunctions.h \
    ../Eigen/src/Core/arch/ZVector/PacketMath.h \
    ../Eigen/src/Core/functors/AssignmentFunctors.h \
    ../Eigen/src/Core/functors/BinaryFunctors.h \
    ../Eigen/src/Core/functors/NullaryFunctors.h \
    ../Eigen/src/Core/functors/StlFunctors.h \
    ../Eigen/src/Core/functors/TernaryFunctors.h \
    ../Eigen/src/Core/functors/UnaryFunctors.h \
    ../Eigen/src/Core/products/GeneralBlockPanelKernel.h \
    ../Eigen/src/Core/products/GeneralMatrixMatrix.h \
    ../Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h \
    ../Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h \
    ../Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
    ../Eigen/src/Core/products/GeneralMatrixVector.h \
    ../Eigen/src/Core/products/GeneralMatrixVector_BLAS.h \
    ../Eigen/src/Core/products/Parallelizer.h \
    ../Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
    ../Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
    ../Eigen/src/Core/products/SelfadjointMatrixVector.h \
    ../Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h \
    ../Eigen/src/Core/products/SelfadjointProduct.h \
    ../Eigen/src/Core/products/SelfadjointRank2Update.h \
    ../Eigen/src/Core/products/TriangularMatrixMatrix.h \
    ../Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h \
    ../Eigen/src/Core/products/TriangularMatrixVector.h \
    ../Eigen/src/Core/products/TriangularMatrixVector_BLAS.h \
    ../Eigen/src/Core/products/TriangularSolverMatrix.h \
    ../Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h \
    ../Eigen/src/Core/products/TriangularSolverVector.h \
    ../Eigen/src/Core/util/BlasUtil.h \
    ../Eigen/src/Core/util/Constants.h \
    ../Eigen/src/Core/util/DisableStupidWarnings.h \
    ../Eigen/src/Core/util/ForwardDeclarations.h \
    ../Eigen/src/Core/util/Macros.h \
    ../Eigen/src/Core/util/Memory.h \
    ../Eigen/src/Core/util/Meta.h \
    ../Eigen/src/Core/util/MKL_support.h \
    ../Eigen/src/Core/util/NonMPL2.h \
    ../Eigen/src/Core/util/ReenableStupidWarnings.h \
    ../Eigen/src/Core/util/StaticAssert.h \
    ../Eigen/src/Core/util/XprHelper.h \
    ../Eigen/src/Core/Array.h \
    ../Eigen/src/Core/ArrayBase.h \
    ../Eigen/src/Core/ArrayWrapper.h \
    ../Eigen/src/Core/Assign.h \
    ../Eigen/src/Core/Assign_MKL.h \
    ../Eigen/src/Core/AssignEvaluator.h \
    ../Eigen/src/Core/BandMatrix.h \
    ../Eigen/src/Core/Block.h \
    ../Eigen/src/Core/BooleanRedux.h \
    ../Eigen/src/Core/CommaInitializer.h \
    ../Eigen/src/Core/ConditionEstimator.h \
    ../Eigen/src/Core/CoreEvaluators.h \
    ../Eigen/src/Core/CoreIterators.h \
    ../Eigen/src/Core/CwiseBinaryOp.h \
    ../Eigen/src/Core/CwiseNullaryOp.h \
    ../Eigen/src/Core/CwiseTernaryOp.h \
    ../Eigen/src/Core/CwiseUnaryOp.h \
    ../Eigen/src/Core/CwiseUnaryView.h \
    ../Eigen/src/Core/DenseBase.h \
    ../Eigen/src/Core/DenseCoeffsBase.h \
    ../Eigen/src/Core/DenseStorage.h \
    ../Eigen/src/Core/Diagonal.h \
    ../Eigen/src/Core/DiagonalMatrix.h \
    ../Eigen/src/Core/DiagonalProduct.h \
    ../Eigen/src/Core/Dot.h \
    ../Eigen/src/Core/EigenBase.h \
    ../Eigen/src/Core/ForceAlignedAccess.h \
    ../Eigen/src/Core/Fuzzy.h \
    ../Eigen/src/Core/GeneralProduct.h \
    ../Eigen/src/Core/GenericPacketMath.h \
    ../Eigen/src/Core/GlobalFunctions.h \
    ../Eigen/src/Core/Inverse.h \
    ../Eigen/src/Core/IO.h \
    ../Eigen/src/Core/Map.h \
    ../Eigen/src/Core/MapBase.h \
    ../Eigen/src/Core/MathFunctions.h \
    ../Eigen/src/Core/MathFunctionsImpl.h \
    ../Eigen/src/Core/Matrix.h \
    ../Eigen/src/Core/MatrixBase.h \
    ../Eigen/src/Core/NestByValue.h \
    ../Eigen/src/Core/NoAlias.h \
    ../Eigen/src/Core/NumTraits.h \
    ../Eigen/src/Core/PermutationMatrix.h \
    ../Eigen/src/Core/PlainObjectBase.h \
    ../Eigen/src/Core/Product.h \
    ../Eigen/src/Core/ProductEvaluators.h \
    ../Eigen/src/Core/Random.h \
    ../Eigen/src/Core/Redux.h \
    ../Eigen/src/Core/Ref.h \
    ../Eigen/src/Core/Replicate.h \
    ../Eigen/src/Core/ReturnByValue.h \
    ../Eigen/src/Core/Reverse.h \
    ../Eigen/src/Core/Select.h \
    ../Eigen/src/Core/SelfAdjointView.h \
    ../Eigen/src/Core/SelfCwiseBinaryOp.h \
    ../Eigen/src/Core/Solve.h \
    ../Eigen/src/Core/SolverBase.h \
    ../Eigen/src/Core/SolveTriangular.h \
    ../Eigen/src/Core/StableNorm.h \
    ../Eigen/src/Core/Stride.h \
    ../Eigen/src/Core/Swap.h \
    ../Eigen/src/Core/Transpose.h \
    ../Eigen/src/Core/Transpositions.h \
    ../Eigen/src/Core/TriangularMatrix.h \
    ../Eigen/src/Core/VectorBlock.h \
    ../Eigen/src/Core/VectorwiseOp.h \
    ../Eigen/src/Core/Visitor.h \
    ../Eigen/src/Eigenvalues/ComplexEigenSolver.h \
    ../Eigen/src/Eigenvalues/ComplexSchur.h \
    ../Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h \
    ../Eigen/src/Eigenvalues/EigenSolver.h \
    ../Eigen/src/Eigenvalues/GeneralizedEigenSolver.h \
    ../Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h \
    ../Eigen/src/Eigenvalues/HessenbergDecomposition.h \
    ../Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h \
    ../Eigen/src/Eigenvalues/RealQZ.h \
    ../Eigen/src/Eigenvalues/RealSchur.h \
    ../Eigen/src/Eigenvalues/RealSchur_LAPACKE.h \
    ../Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h \
    ../Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h \
    ../Eigen/src/Eigenvalues/Tridiagonalization.h \
    ../Eigen/src/Geometry/arch/Geometry_SSE.h \
    ../Eigen/src/Geometry/AlignedBox.h \
    ../Eigen/src/Geometry/AngleAxis.h \
    ../Eigen/src/Geometry/EulerAngles.h \
    ../Eigen/src/Geometry/Homogeneous.h \
    ../Eigen/src/Geometry/Hyperplane.h \
    ../Eigen/src/Geometry/OrthoMethods.h \
    ../Eigen/src/Geometry/ParametrizedLine.h \
    ../Eigen/src/Geometry/Quaternion.h \
    ../Eigen/src/Geometry/Rotation2D.h \
    ../Eigen/src/Geometry/RotationBase.h \
    ../Eigen/src/Geometry/Scaling.h \
    ../Eigen/src/Geometry/Transform.h \
    ../Eigen/src/Geometry/Translation.h \
    ../Eigen/src/Geometry/Umeyama.h \
    ../Eigen/src/Householder/BlockHouseholder.h \
    ../Eigen/src/Householder/Householder.h \
    ../Eigen/src/Householder/HouseholderSequence.h \
    ../Eigen/src/Jacobi/Jacobi.h \
    ../Eigen/src/LU/arch/Inverse_SSE.h \
    ../Eigen/src/LU/Determinant.h \
    ../Eigen/src/LU/FullPivLU.h \
    ../Eigen/src/LU/InverseImpl.h \
    ../Eigen/src/LU/PartialPivLU.h \
    ../Eigen/src/LU/PartialPivLU_LAPACKE.h \
    ../Eigen/src/misc/blas.h \
    ../Eigen/src/misc/Image.h \
    ../Eigen/src/misc/Kernel.h \
    ../Eigen/src/misc/lapack.h \
    ../Eigen/src/misc/lapacke.h \
    ../Eigen/src/misc/lapacke_mangling.h \
    ../Eigen/src/misc/RealSvd2x2.h \
    ../Eigen/src/plugins/ArrayCwiseBinaryOps.h \
    ../Eigen/src/plugins/ArrayCwiseUnaryOps.h \
    ../Eigen/src/plugins/BlockMethods.h \
    ../Eigen/src/plugins/CommonCwiseBinaryOps.h \
    ../Eigen/src/plugins/CommonCwiseUnaryOps.h \
    ../Eigen/src/plugins/MatrixCwiseBinaryOps.h \
    ../Eigen/src/plugins/MatrixCwiseUnaryOps.h \
    ../Eigen/src/QR/ColPivHouseholderQR.h \
    ../Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h \
    ../Eigen/src/QR/CompleteOrthogonalDecomposition.h \
    ../Eigen/src/QR/FullPivHouseholderQR.h \
    ../Eigen/src/QR/HouseholderQR.h \
    ../Eigen/src/QR/HouseholderQR_LAPACKE.h \
    ../Eigen/src/SVD/BDCSVD.h \
    ../Eigen/src/SVD/JacobiSVD.h \
    ../Eigen/src/SVD/JacobiSVD_LAPACKE.h \
    ../Eigen/src/SVD/SVDBase.h \
    ../Eigen/src/SVD/UpperBidiagonalization.h \
    ../Parser/MemoryMapped/MemoryMapped.h \
    ../Variant.h \
    ../Output/OutputHandler.h \
    ../Request.h \
    ../Parser/File.h \
    ../vikNGS.h \
    src/windows/MainWindow.h \
    src/windows/PlotWindow.h \
    src/log/TypeConverter.h \
    src/AsyncJob.h \
    ../SampleInfo.h \
    ../Parser/Parser.h \
    ../Parser/Filter.h \
    ../Math/Math.h \
    ../Interval.h \
    ../Test/Test.h \
    ../Test/TestObject.h \
    src/windows/Chromosome.h \
    ../Log.h \
    ../Math/CompQuadForm.h \
    src/windows/TableDisplayWindow.h \
    src/simulation/Simulation.h \
    src/log/qlog.h \
    src/widgets/qzoombar.h \
    src/widgets/qcustomplot.h \
    src/windows/SimPlotWindow.h \
    src/widgets/qcpdocumentobject.h \
    ../Test/ScoreTestFunctions.h \
    ../Test/Group.h \
    ../Test/Genotype.h \
    ../Test/Phenotype.h \
    src/windows/VikngsUiConstants.h \
    ../Enum/Test.h \
    ../Enum/CollapseType.h \
    ../Enum/Depth.h \
    ../Enum/Family.h \
    ../Enum/Filter.h \
    ../Enum/GenotypeSource.h \
    ../Enum/Statistic.h \
    ../Enum/Test.h \
    ../Enum/Variance.h


FORMS += \
    ui/mainwindow.ui \
    ui/plotwindow.ui \
    ui/simplotwindow.ui \
    ui/tabledisplaywindow.ui


