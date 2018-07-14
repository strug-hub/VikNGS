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
    src/vikngs.cpp \
    src/windows/maintab.cpp \
    src/windows/plotwindow.cpp \
    src/widgets/qzoombar.cpp \
    src/widgets/qcustomplot.cpp \
    src/log/log.cpp \
    src/windows/plotwindowmouse.cpp \
    src/windows/simplotwindow.cpp \
    src/simulation/simulation.cpp \
    src/simulation/simulationhelper.cpp \
    ../src/Test/CommonTest.cpp \
    ../src/Test/RareTest.cpp \
    ../src/Test/RareTestObject.cpp \
    ../src/RequestBuilder.cpp \
    ../src/Math/StatisticsHelper.cpp \
    ../src/Math/CompQuadForm.cpp \
    ../src/Math/VectorHelper.cpp \
    ../src/Parser/MemoryMapped/MemoryMapped.cpp \
    ../src/Parser/VCF/VCFParserUtils.cpp \
    ../src/Parser/Filter/VariantFilter.cpp \
    ../src/Parser/Sample/SampleParserUtils.cpp \
    ../src/Parser/Sample/SampleParser.cpp \
    ../src/Parser/BED/BEDParserUtils.cpp \
    ../src/Parser/BED/BEDParser.cpp \
    ../src/Parser/InputParser.cpp \
    ../src/Parser/ParserTools.cpp \
    ../src/Test/CommonTestObject.cpp \
    src/windows/simulationtab.cpp \
    ../src/Test/RareTestCollapseObject.cpp \
    src/global.cpp \
    src/windows/qsimulationtab.cpp

HEADERS += \
    src/widgets/qcustomplot.h \
    src/windows/mainwindow.h \
    src/runner.h \
    src/log/qlog.h \
    src/log/typeconverter.h \
    src/windows/plotwindow.h \
    src/widgets/qzoombar.h \
    src/windows/simplotwindow.h \
    src/simulation/simulation.h \
    ../src/RVS.h \
    ../src/Variant.h \
    ../src/Eigen/src/Cholesky/LDLT.h \
    ../src/Eigen/src/Cholesky/LLT.h \
    ../src/Eigen/src/Cholesky/LLT_LAPACKE.h \
    ../src/Eigen/src/Core/arch/AltiVec/Complex.h \
    ../src/Eigen/src/Core/arch/AltiVec/MathFunctions.h \
    ../src/Eigen/src/Core/arch/AltiVec/PacketMath.h \
    ../src/Eigen/src/Core/arch/AVX/Complex.h \
    ../src/Eigen/src/Core/arch/AVX/MathFunctions.h \
    ../src/Eigen/src/Core/arch/AVX/PacketMath.h \
    ../src/Eigen/src/Core/arch/AVX/TypeCasting.h \
    ../src/Eigen/src/Core/arch/AVX512/MathFunctions.h \
    ../src/Eigen/src/Core/arch/AVX512/PacketMath.h \
    ../src/Eigen/src/Core/arch/CUDA/Complex.h \
    ../src/Eigen/src/Core/arch/CUDA/Half.h \
    ../src/Eigen/src/Core/arch/CUDA/MathFunctions.h \
    ../src/Eigen/src/Core/arch/CUDA/PacketMath.h \
    ../src/Eigen/src/Core/arch/CUDA/PacketMathHalf.h \
    ../src/Eigen/src/Core/arch/CUDA/TypeCasting.h \
    ../src/Eigen/src/Core/arch/Default/Settings.h \
    ../src/Eigen/src/Core/arch/NEON/Complex.h \
    ../src/Eigen/src/Core/arch/NEON/MathFunctions.h \
    ../src/Eigen/src/Core/arch/NEON/PacketMath.h \
    ../src/Eigen/src/Core/arch/SSE/Complex.h \
    ../src/Eigen/src/Core/arch/SSE/MathFunctions.h \
    ../src/Eigen/src/Core/arch/SSE/PacketMath.h \
    ../src/Eigen/src/Core/arch/SSE/TypeCasting.h \
    ../src/Eigen/src/Core/arch/ZVector/Complex.h \
    ../src/Eigen/src/Core/arch/ZVector/MathFunctions.h \
    ../src/Eigen/src/Core/arch/ZVector/PacketMath.h \
    ../src/Eigen/src/Core/functors/AssignmentFunctors.h \
    ../src/Eigen/src/Core/functors/BinaryFunctors.h \
    ../src/Eigen/src/Core/functors/NullaryFunctors.h \
    ../src/Eigen/src/Core/functors/StlFunctors.h \
    ../src/Eigen/src/Core/functors/TernaryFunctors.h \
    ../src/Eigen/src/Core/functors/UnaryFunctors.h \
    ../src/Eigen/src/Core/products/GeneralBlockPanelKernel.h \
    ../src/Eigen/src/Core/products/GeneralMatrixMatrix.h \
    ../src/Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h \
    ../src/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h \
    ../src/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
    ../src/Eigen/src/Core/products/GeneralMatrixVector.h \
    ../src/Eigen/src/Core/products/GeneralMatrixVector_BLAS.h \
    ../src/Eigen/src/Core/products/Parallelizer.h \
    ../src/Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
    ../src/Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
    ../src/Eigen/src/Core/products/SelfadjointMatrixVector.h \
    ../src/Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h \
    ../src/Eigen/src/Core/products/SelfadjointProduct.h \
    ../src/Eigen/src/Core/products/SelfadjointRank2Update.h \
    ../src/Eigen/src/Core/products/TriangularMatrixMatrix.h \
    ../src/Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h \
    ../src/Eigen/src/Core/products/TriangularMatrixVector.h \
    ../src/Eigen/src/Core/products/TriangularMatrixVector_BLAS.h \
    ../src/Eigen/src/Core/products/TriangularSolverMatrix.h \
    ../src/Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h \
    ../src/Eigen/src/Core/products/TriangularSolverVector.h \
    ../src/Eigen/src/Core/util/BlasUtil.h \
    ../src/Eigen/src/Core/util/Constants.h \
    ../src/Eigen/src/Core/util/DisableStupidWarnings.h \
    ../src/Eigen/src/Core/util/ForwardDeclarations.h \
    ../src/Eigen/src/Core/util/Macros.h \
    ../src/Eigen/src/Core/util/Memory.h \
    ../src/Eigen/src/Core/util/Meta.h \
    ../src/Eigen/src/Core/util/MKL_support.h \
    ../src/Eigen/src/Core/util/NonMPL2.h \
    ../src/Eigen/src/Core/util/ReenableStupidWarnings.h \
    ../src/Eigen/src/Core/util/StaticAssert.h \
    ../src/Eigen/src/Core/util/XprHelper.h \
    ../src/Eigen/src/Core/Array.h \
    ../src/Eigen/src/Core/ArrayBase.h \
    ../src/Eigen/src/Core/ArrayWrapper.h \
    ../src/Eigen/src/Core/Assign.h \
    ../src/Eigen/src/Core/Assign_MKL.h \
    ../src/Eigen/src/Core/AssignEvaluator.h \
    ../src/Eigen/src/Core/BandMatrix.h \
    ../src/Eigen/src/Core/Block.h \
    ../src/Eigen/src/Core/BooleanRedux.h \
    ../src/Eigen/src/Core/CommaInitializer.h \
    ../src/Eigen/src/Core/ConditionEstimator.h \
    ../src/Eigen/src/Core/CoreEvaluators.h \
    ../src/Eigen/src/Core/CoreIterators.h \
    ../src/Eigen/src/Core/CwiseBinaryOp.h \
    ../src/Eigen/src/Core/CwiseNullaryOp.h \
    ../src/Eigen/src/Core/CwiseTernaryOp.h \
    ../src/Eigen/src/Core/CwiseUnaryOp.h \
    ../src/Eigen/src/Core/CwiseUnaryView.h \
    ../src/Eigen/src/Core/DenseBase.h \
    ../src/Eigen/src/Core/DenseCoeffsBase.h \
    ../src/Eigen/src/Core/DenseStorage.h \
    ../src/Eigen/src/Core/Diagonal.h \
    ../src/Eigen/src/Core/DiagonalMatrix.h \
    ../src/Eigen/src/Core/DiagonalProduct.h \
    ../src/Eigen/src/Core/Dot.h \
    ../src/Eigen/src/Core/EigenBase.h \
    ../src/Eigen/src/Core/ForceAlignedAccess.h \
    ../src/Eigen/src/Core/Fuzzy.h \
    ../src/Eigen/src/Core/GeneralProduct.h \
    ../src/Eigen/src/Core/GenericPacketMath.h \
    ../src/Eigen/src/Core/GlobalFunctions.h \
    ../src/Eigen/src/Core/Inverse.h \
    ../src/Eigen/src/Core/IO.h \
    ../src/Eigen/src/Core/Map.h \
    ../src/Eigen/src/Core/MapBase.h \
    ../src/Eigen/src/Core/MathFunctions.h \
    ../src/Eigen/src/Core/MathFunctionsImpl.h \
    ../src/Eigen/src/Core/Matrix.h \
    ../src/Eigen/src/Core/MatrixBase.h \
    ../src/Eigen/src/Core/NestByValue.h \
    ../src/Eigen/src/Core/NoAlias.h \
    ../src/Eigen/src/Core/NumTraits.h \
    ../src/Eigen/src/Core/PermutationMatrix.h \
    ../src/Eigen/src/Core/PlainObjectBase.h \
    ../src/Eigen/src/Core/Product.h \
    ../src/Eigen/src/Core/ProductEvaluators.h \
    ../src/Eigen/src/Core/Random.h \
    ../src/Eigen/src/Core/Redux.h \
    ../src/Eigen/src/Core/Ref.h \
    ../src/Eigen/src/Core/Replicate.h \
    ../src/Eigen/src/Core/ReturnByValue.h \
    ../src/Eigen/src/Core/Reverse.h \
    ../src/Eigen/src/Core/Select.h \
    ../src/Eigen/src/Core/SelfAdjointView.h \
    ../src/Eigen/src/Core/SelfCwiseBinaryOp.h \
    ../src/Eigen/src/Core/Solve.h \
    ../src/Eigen/src/Core/SolverBase.h \
    ../src/Eigen/src/Core/SolveTriangular.h \
    ../src/Eigen/src/Core/StableNorm.h \
    ../src/Eigen/src/Core/Stride.h \
    ../src/Eigen/src/Core/Swap.h \
    ../src/Eigen/src/Core/Transpose.h \
    ../src/Eigen/src/Core/Transpositions.h \
    ../src/Eigen/src/Core/TriangularMatrix.h \
    ../src/Eigen/src/Core/VectorBlock.h \
    ../src/Eigen/src/Core/VectorwiseOp.h \
    ../src/Eigen/src/Core/Visitor.h \
    ../src/Eigen/src/Eigenvalues/ComplexEigenSolver.h \
    ../src/Eigen/src/Eigenvalues/ComplexSchur.h \
    ../src/Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h \
    ../src/Eigen/src/Eigenvalues/EigenSolver.h \
    ../src/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h \
    ../src/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h \
    ../src/Eigen/src/Eigenvalues/HessenbergDecomposition.h \
    ../src/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h \
    ../src/Eigen/src/Eigenvalues/RealQZ.h \
    ../src/Eigen/src/Eigenvalues/RealSchur.h \
    ../src/Eigen/src/Eigenvalues/RealSchur_LAPACKE.h \
    ../src/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h \
    ../src/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h \
    ../src/Eigen/src/Eigenvalues/Tridiagonalization.h \
    ../src/Eigen/src/Geometry/arch/Geometry_SSE.h \
    ../src/Eigen/src/Geometry/AlignedBox.h \
    ../src/Eigen/src/Geometry/AngleAxis.h \
    ../src/Eigen/src/Geometry/EulerAngles.h \
    ../src/Eigen/src/Geometry/Homogeneous.h \
    ../src/Eigen/src/Geometry/Hyperplane.h \
    ../src/Eigen/src/Geometry/OrthoMethods.h \
    ../src/Eigen/src/Geometry/ParametrizedLine.h \
    ../src/Eigen/src/Geometry/Quaternion.h \
    ../src/Eigen/src/Geometry/Rotation2D.h \
    ../src/Eigen/src/Geometry/RotationBase.h \
    ../src/Eigen/src/Geometry/Scaling.h \
    ../src/Eigen/src/Geometry/Transform.h \
    ../src/Eigen/src/Geometry/Translation.h \
    ../src/Eigen/src/Geometry/Umeyama.h \
    ../src/Eigen/src/Householder/BlockHouseholder.h \
    ../src/Eigen/src/Householder/Householder.h \
    ../src/Eigen/src/Householder/HouseholderSequence.h \
    ../src/Eigen/src/Jacobi/Jacobi.h \
    ../src/Eigen/src/LU/arch/Inverse_SSE.h \
    ../src/Eigen/src/LU/Determinant.h \
    ../src/Eigen/src/LU/FullPivLU.h \
    ../src/Eigen/src/LU/InverseImpl.h \
    ../src/Eigen/src/LU/PartialPivLU.h \
    ../src/Eigen/src/LU/PartialPivLU_LAPACKE.h \
    ../src/Eigen/src/misc/blas.h \
    ../src/Eigen/src/misc/Image.h \
    ../src/Eigen/src/misc/Kernel.h \
    ../src/Eigen/src/misc/lapack.h \
    ../src/Eigen/src/misc/lapacke.h \
    ../src/Eigen/src/misc/lapacke_mangling.h \
    ../src/Eigen/src/misc/RealSvd2x2.h \
    ../src/Eigen/src/plugins/ArrayCwiseBinaryOps.h \
    ../src/Eigen/src/plugins/ArrayCwiseUnaryOps.h \
    ../src/Eigen/src/plugins/BlockMethods.h \
    ../src/Eigen/src/plugins/CommonCwiseBinaryOps.h \
    ../src/Eigen/src/plugins/CommonCwiseUnaryOps.h \
    ../src/Eigen/src/plugins/MatrixCwiseBinaryOps.h \
    ../src/Eigen/src/plugins/MatrixCwiseUnaryOps.h \
    ../src/Eigen/src/QR/ColPivHouseholderQR.h \
    ../src/Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h \
    ../src/Eigen/src/QR/CompleteOrthogonalDecomposition.h \
    ../src/Eigen/src/QR/FullPivHouseholderQR.h \
    ../src/Eigen/src/QR/HouseholderQR.h \
    ../src/Eigen/src/QR/HouseholderQR_LAPACKE.h \
    ../src/Eigen/src/SVD/BDCSVD.h \
    ../src/Eigen/src/SVD/JacobiSVD.h \
    ../src/Eigen/src/SVD/JacobiSVD_LAPACKE.h \
    ../src/Eigen/src/SVD/SVDBase.h \
    ../src/Eigen/src/SVD/UpperBidiagonalization.h \
    ../src/Parser/MemoryMapped/MemoryMapped.h \
    ../src/Parser/InputParser.h \
    ../src/Output/OutputHandler.h \
    ../src/Log.h \
    ../src/Request.h \
    ../src/RVS.h \
    ../src/Math/MathHelper.h \
    ../src/Math/StatisticsHelper.h \
    ../src/Math/VectorHelper.h \
    ../src/Parser/VCF/VCFParserUtils.h \
    ../src/Parser/Sample/SampleParserUtils.h \
    ../src/Parser/BED/BEDParserUtils.h \
    ../src/Parser/Filter/VariantFilterUtils.h \
    ../src/Parser/File.h \
    ../src/Test/CommonTestObject.h \
    ../src/Test/RareTestObject.h \
    ../src/Test/RareTestCollapseObject.h \
    ../src/Parser/BED/Interval.h \
    src/global.h


FORMS += \
    ui/mainwindow.ui \
    ui/plotwindow.ui \
    ui/simplotwindow.ui

