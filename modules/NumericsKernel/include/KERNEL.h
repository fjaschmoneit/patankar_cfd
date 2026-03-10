#include <unordered_map>
#include <variant>
#include <memory>
#include "KernelTypeDefs.h"

namespace KERNEL {
    enum SolverMethod {
        GaussSeidel, BiCGSTAB, Jacobi, Blaze_automatic
    };


    smatrix newTempBandedSMatrix(std::size_t N, const std::vector<int>& bandIDs, GLOBAL::scalar init);


    // test performance of these
    template<typename MatrixType>
    void fillBand(blaze::Band<MatrixType> band, GLOBAL::scalar value) {
        for(size_t i = 0; i < band.size(); i++)  band[i] = value;
    }

    template<typename MatrixType>
    void fillBand(blaze::Band<MatrixType> band, vector &values) {
        // FJA test if length values matches length band
        for(size_t i = 0; i < band.size(); i++)  band[i] = values[i];
    }

    std::vector<int> getBandIDs(const smatrix &A);

    void solve(const smatrix& A, vector& x, const vector& b, SolverMethod method, const GLOBAL::scalar tolerance=1e-15, const unsigned int maxIter=10000);
    void solve(const dmatrix& A, vector& x, const vector& b, SolverMethod method, const GLOBAL::scalar tolerance=1e-15, const unsigned int maxIter=10000);

    // Variant holding unique_ptrs to vector or matrix
    using ObjectVariantPtr = std::variant<
        std::unique_ptr<vector>,
        std::unique_ptr<dmatrix>,
        std::unique_ptr<smatrix>
    >;

    struct VectorHandle { size_t id; };
    struct MatrixHandle { size_t id; };

    class ObjectRegistry {
    private:
        std::unordered_map<size_t, ObjectVariantPtr> registry_;
        size_t nextID_ = 0;
        bool registryClosed_ = false;

    public:

        // I am closing the reg after initial problem setup.
        // I can therefore savely return references to its objects
        // from associated getter functions
        void closeRegistry() {
            registryClosed_ = true;
            registry_.reserve(registry_.size());  // Prevent reallocation
        }

        // definition must be together with blaze includes.
        // Here, together with blaz's forward declaration,
        // I can only declare:
        ObjectRegistry();
        ~ObjectRegistry();
        ObjectRegistry(ObjectRegistry&&) noexcept;
        ObjectRegistry& operator=(ObjectRegistry&&) noexcept;

        // deleting copy constructor:
        ObjectRegistry(const ObjectRegistry&) = delete;
        ObjectRegistry& operator=(const ObjectRegistry&) = delete;

        VectorHandle newVector(size_t size, GLOBAL::scalar initialValue = 0.0);
        MatrixHandle newMatrix(size_t rows, size_t cols, bool sparse = false );


        // ============ Getters ===================

        // these object references are valid, because I cannot create or
        // delete new objects in the registry after it is closed.

        KERNEL::vector& getVectorRef(VectorHandle handle);

        KERNEL::dmatrix& getDenseMatrixRef(MatrixHandle handle);
        KERNEL::smatrix& getSparseMatrixRef(MatrixHandle handle);

    };
}