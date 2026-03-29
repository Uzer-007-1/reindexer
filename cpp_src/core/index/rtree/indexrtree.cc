#include "indexrtree.h"

#include "core/rdxcontext.h"
#include <unordered_set>
#include "greenesplitter.h"
#include "linearsplitter.h"
#include "quadraticsplitter.h"
#include "rstarsplitter.h"

namespace reindexer {

namespace {
bool isValidLongitude(double x) noexcept { return x >= -180.0 && x <= 180.0; }
bool isValidLatitude(double y) noexcept { return y >= -90.0 && y <= 90.0; }
double metersToDegrees(double distanceMeters) noexcept {
	constexpr double kMetersPerDegreeAtEquator = 111319.49079327358;
	return distanceMeters / kMetersPerDegreeAtEquator;
}
void validateGeoPoint(Point p, std::string_view context) {
	if (!isValidLongitude(p.X()) || !isValidLatitude(p.Y())) {
		throw Error(errParams, "{}: GIS point is out of WGS84 range. Expected [lon, lat], got [{}, {}]", context, p.X(), p.Y());
	}
}
void validateGeoIndexPoint(Point p, std::string_view context) {
	if (!isValidLongitude(p.X()) || !isValidLatitude(p.Y())) {
		throw Error(errQueryExec, "{}: GIS index contains point out of WGS84 range. Expected [lon, lat], got [{}, {}]", context, p.X(),
					p.Y());
	}
}
}  // namespace

template <typename KeyEntryT, template <typename, typename, typename, typename, size_t, size_t> class Splitter, size_t MaxEntries,
		  size_t MinEntries>
SelectKeyResults IndexRTree<KeyEntryT, Splitter, MaxEntries, MinEntries>::SelectKey(const VariantArray& keys, CondType condition,
																					SortType sortId, const Index::SelectContext& selectCtx,
																					const RdxContext& rdxCtx) {
	const auto indexWard(rdxCtx.BeforeIndexWork());
	// Generic comparator path is planar and should not be used for GIS semantics.
	if (selectCtx.opts.forceComparator && !this->opts_.IsGeo()) {
		return IndexStore<typename Map::key_type>::SelectKey(keys, condition, sortId, selectCtx, rdxCtx);
	}

	SelectKeyResult res;

	if (condition != CondDWithin) {
		throw Error(errQueryExec, "Only DWithin condition is available for RTree index");
	}
	if (keys.size() != 2) {
		throw Error(errQueryExec, "DWithin condition expects two arguments");
	}
	Point point;
	double distance;
	if (keys[0].Type().Is<KeyValueType::Tuple>()) {
		point = keys[0].As<Point>();
		distance = keys[1].As<double>();
	} else {
		point = keys[1].As<Point>();
		distance = keys[0].As<double>();
	}
	if (this->opts_.IsGeo()) {
		validateGeoPoint(point, "Select DWithin");
		if (distance < 0.0) {
			throw Error(errParams, "Select DWithin: GIS distance can not be negative");
		}
	}
	class [[nodiscard]] Visitor : public Map::Visitor {
	public:
		Visitor(SortType sId, unsigned distinct, unsigned iCountInNs, SelectKeyResult& r)
			: sortId_{sId}, itemsCountInNs_{distinct ? 0u : iCountInNs}, res_{r} {}
		bool operator()(const typename Map::value_type& v) override {
			idsCount_ += v.second.Unsorted().size();
			res_.emplace_back(v.second, sortId_);
			return ScanWin();
		}
		bool ScanWin() const noexcept {
			return itemsCountInNs_ && res_.size() > 1u && (100u * idsCount_ / itemsCountInNs_ > maxSelectivityPercentForIdset());
		}

	private:
		SortType sortId_;
		unsigned itemsCountInNs_;
		SelectKeyResult& res_;
		size_t idsCount_ = 0;
	} visitor{sortId, selectCtx.opts.distinct, selectCtx.opts.itemsCountInNamespace, res};
	if (this->opts_.IsGeo()) {
		class [[nodiscard]] GeoVisitor : public Map::Visitor {
		public:
			GeoVisitor(SortType sId, unsigned distinct, unsigned iCountInNs, SelectKeyResult& r, Point c, double dMeters)
				: sortId_{sId}, itemsCountInNs_{distinct ? 0u : iCountInNs}, res_{r}, center_{c}, distanceMeters_{dMeters} {}
			bool operator()(const typename Map::value_type& v) override {
				validateGeoIndexPoint(v.first, "Select DWithin");
				if (!visitedPoints_.emplace(v.first).second) {
					return false;
				}
				if (!GeoDWithin(v.first, center_, distanceMeters_)) {
					return false;
				}
				idsCount_ += v.second.Unsorted().size();
				res_.emplace_back(v.second, sortId_);
				return ScanWin();
			}
			bool ScanWin() const noexcept {
				return itemsCountInNs_ && res_.size() > 1u && (100u * idsCount_ / itemsCountInNs_ > maxSelectivityPercentForIdset());
			}

		private:
			SortType sortId_;
			unsigned itemsCountInNs_;
			SelectKeyResult& res_;
			Point center_;
			double distanceMeters_;
			std::unordered_set<Point, std::hash<Point>, point_strict_equal> visitedPoints_;
			size_t idsCount_ = 0;
		} geoVisitor{sortId, selectCtx.opts.distinct, selectCtx.opts.itemsCountInNamespace, res, point, distance};
		const double degRadius = metersToDegrees(distance);
		this->idx_map.DWithin(point, degRadius, geoVisitor);
		// Handle antimeridian wrap-around for GIS longitudes.
		const double minLon = point.X() - degRadius;
		const double maxLon = point.X() + degRadius;
		if (minLon < -180.0) {
			this->idx_map.DWithin(Point{point.X() + 360.0, point.Y()}, degRadius, geoVisitor);
		}
		if (maxLon > 180.0) {
			this->idx_map.DWithin(Point{point.X() - 360.0, point.Y()}, degRadius, geoVisitor);
		}
		// NOTE: Fallback to generic comparator is intentionally disabled for GIS mode,
		// because generic DWithin comparator works in planar coordinates and may produce
		// incorrect results for geodesic (meters-based) semantics.
		return SelectKeyResults(std::move(res));
	}
	this->idx_map.DWithin(point, distance, visitor);
	if (visitor.ScanWin()) {
		// fallback to comparator, due to expensive idset
		return IndexStore<typename Map::key_type>::SelectKey(keys, condition, sortId, selectCtx, rdxCtx);
	}
	return SelectKeyResults(std::move(res));
}

template <typename KeyEntryT, template <typename, typename, typename, typename, size_t, size_t> class Splitter, size_t MaxEntries,
		  size_t MinEntries>
void IndexRTree<KeyEntryT, Splitter, MaxEntries, MinEntries>::Upsert(VariantArray& result, const VariantArray& keys, IdType id,
																	 bool& clearCache) {
	if (keys.empty() || keys.IsNullValue()) {
		std::ignore = Upsert(Variant{}, id, clearCache);
		return;
	}
	if (keys.size() != 2) [[unlikely]] {
		[[maybe_unused]] Point p{keys};
	}
	const Point point{keys[0].IsNullValue() ? 0.0 : keys[0].As<double>(), keys[1].IsNullValue() ? 0.0 : keys[1].As<double>()};
	if (this->opts_.IsGeo()) {
		validateGeoPoint(point, "Upsert");
	}
	typename Map::iterator keyIt = this->idx_map.find(point);
	if (keyIt == this->idx_map.end()) {
		keyIt = this->idx_map.insert_without_test({point, typename Map::mapped_type()});
	} else {
		this->delMemStat(keyIt);
	}

	if (keyIt->second.Unsorted().Add(id, this->opts_.IsPK() ? IdSetEditMode::Ordered : IdSetEditMode::Auto, this->sortedIdxCount_)) {
		this->isBuilt_ = false;
		// reset cache
		this->cache_.ResetImpl();
		clearCache = true;
	}
	this->tracker_.markUpdated(this->idx_map, keyIt);

	this->addMemStat(keyIt);

	result = VariantArray(keyIt->first);
}

template <typename KeyEntryT, template <typename, typename, typename, typename, size_t, size_t> class Splitter, size_t MaxEntries,
		  size_t MinEntries>
void IndexRTree<KeyEntryT, Splitter, MaxEntries, MinEntries>::Delete(const VariantArray& keys, IdType id, MustExist mustExist,
																	 StringsHolder& strHolder, bool& clearCache) {
	if (keys.empty() || keys.IsNullValue()) {
		return Delete(Variant{}, id, mustExist, strHolder, clearCache);
	}
	int delcnt = 0;
	const Point point = static_cast<Point>(keys);
	typename Map::iterator keyIt = this->idx_map.find(point);
	if (keyIt == this->idx_map.end()) {
		return;
	}
	this->cache_.ResetImpl();
	clearCache = true;
	this->isBuilt_ = false;

	this->delMemStat(keyIt);
	delcnt = keyIt->second.Unsorted().Erase(id);
	(void)delcnt;
	assertf(!mustExist || this->Opts().IsSparse() || delcnt, "Delete non-existent id from index '{}' id={}, key={} ({})", this->name_, id,
			Variant(keys).template As<std::string>(this->payloadType_, this->Fields()),
			Variant(keyIt->first).As<std::string>(this->payloadType_, this->Fields()));

	if (keyIt->second.Unsorted().IsEmpty()) {
		this->tracker_.markDeleted(keyIt);
		this->idx_map.template erase<void>(keyIt);
	} else {
		this->addMemStat(keyIt);
		this->tracker_.markUpdated(this->idx_map, keyIt);
	}
}

std::unique_ptr<Index> IndexRTree_New(const IndexDef& idef, PayloadType&& payloadType, FieldsSet&& fields,
									  const NamespaceCacheConfigData& cacheCfg) {
	assertrx_throw(!idef.Opts().IsPK());
	switch (idef.Opts().RTreeType()) {
		case IndexOpts::Linear:
			if (idef.Opts().IsDense()) {
				return std::make_unique<IndexRTree<Index::KeyEntryPlain, LinearSplitter, 32, 4>>(idef, std::move(payloadType),
																								 std::move(fields), cacheCfg);
			} else {
				return std::make_unique<IndexRTree<Index::KeyEntry, LinearSplitter, 32, 4>>(idef, std::move(payloadType), std::move(fields),
																							cacheCfg);
			}
		case IndexOpts::Quadratic:
			if (idef.Opts().IsDense()) {
				return std::make_unique<IndexRTree<Index::KeyEntryPlain, QuadraticSplitter, 32, 4>>(idef, std::move(payloadType),
																									std::move(fields), cacheCfg);
			}
			return std::make_unique<IndexRTree<Index::KeyEntry, QuadraticSplitter, 32, 4>>(idef, std::move(payloadType), std::move(fields),
																						   cacheCfg);
		case IndexOpts::Greene:
			if (idef.Opts().IsDense()) {
				return std::make_unique<IndexRTree<Index::KeyEntryPlain, GreeneSplitter, 16, 4>>(idef, std::move(payloadType),
																								 std::move(fields), cacheCfg);
			} else {
				return std::make_unique<IndexRTree<Index::KeyEntry, GreeneSplitter, 16, 4>>(idef, std::move(payloadType), std::move(fields),
																							cacheCfg);
			}
		case IndexOpts::RStar:
			if (idef.Opts().IsDense()) {
				return std::make_unique<IndexRTree<Index::KeyEntryPlain, RStarSplitter, 32, 4>>(idef, std::move(payloadType),
																								std::move(fields), cacheCfg);
			} else {
				return std::make_unique<IndexRTree<Index::KeyEntry, RStarSplitter, 32, 4>>(idef, std::move(payloadType), std::move(fields),
																						   cacheCfg);
			}
	}

	assertrx(0);
	std::abort();
}

template class IndexRTree<Index::KeyEntry, LinearSplitter, 32, 4>;
template class IndexRTree<Index::KeyEntryPlain, LinearSplitter, 32, 4>;
template class IndexRTree<Index::KeyEntry, QuadraticSplitter, 32, 4>;
template class IndexRTree<Index::KeyEntryPlain, QuadraticSplitter, 32, 4>;
template class IndexRTree<Index::KeyEntry, GreeneSplitter, 16, 4>;
template class IndexRTree<Index::KeyEntryPlain, GreeneSplitter, 16, 4>;
template class IndexRTree<Index::KeyEntry, RStarSplitter, 32, 4>;
template class IndexRTree<Index::KeyEntryPlain, RStarSplitter, 32, 4>;

}  // namespace reindexer
