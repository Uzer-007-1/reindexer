package reindexer

import (
	"math"
	"reflect"
	"testing"
)

func TestValidateGeoPoint(t *testing.T) {
	if err := ValidateGeoPoint(Point{0, 0}); err != nil {
		t.Fatalf("expected valid point, got error: %v", err)
	}
	if err := ValidateGeoPoint(Point{-181, 0}); err == nil {
		t.Fatal("expected longitude validation error")
	}
	if err := ValidateGeoPoint(Point{0, 91}); err == nil {
		t.Fatal("expected latitude validation error")
	}
}

func TestGeoDistanceMeters(t *testing.T) {
	london := Point{-0.1278, 51.5074}
	paris := Point{2.3522, 48.8566}
	dist := GeoDistanceMeters(london, paris)
	if dist < 340000 || dist > 350000 {
		t.Fatalf("unexpected london-paris distance: %v", dist)
	}
}

func TestGeoDWithin(t *testing.T) {
	p1 := Point{37.6173, 55.7558} // Moscow
	p2 := Point{30.3351, 59.9343} // Saint Petersburg
	if !GeoDWithin(p1, p2, 700000) {
		t.Fatal("expected points to be within 700km")
	}
	if GeoDWithin(p1, p2, 500000) {
		t.Fatal("expected points to be farther than 500km")
	}
}

func TestGeoDistanceMetersAntipodesIsFinite(t *testing.T) {
	d := GeoDistanceMeters(Point{0, 0}, Point{180, 0})
	if math.IsNaN(d) || math.IsInf(d, 0) {
		t.Fatalf("distance must be finite for antipodes, got: %v", d)
	}
	if d < 20000000 || d > 20160000 {
		t.Fatalf("unexpected antipodal distance: %v", d)
	}
}

func TestDWithinGeoCheckedValidation(t *testing.T) {
	q := &Query{}
	gotQ, err := q.DWithinGeoChecked("pt", Point{181, 0}, 10)
	if err == nil {
		t.Fatal("expected validation error for longitude")
	}
	if gotQ != q {
		t.Fatal("expected original query pointer on validation error")
	}

	gotQ, err = q.DWithinGeoChecked("pt", Point{10, 10}, -1)
	if err == nil {
		t.Fatal("expected validation error for negative distance")
	}
	if gotQ != q {
		t.Fatal("expected original query pointer on validation error")
	}
}

func TestSortStGeoPointDistanceCheckedValidation(t *testing.T) {
	q := &Query{}
	gotQ, err := q.SortStGeoPointDistanceChecked("pt", Point{181, 0}, false)
	if err == nil {
		t.Fatal("expected validation error for longitude")
	}
	if gotQ != q {
		t.Fatal("expected original query pointer on validation error")
	}

	gotQ, err = q.SortStGeoPointDistanceChecked("pt", Point{10, 20}, true)
	if err != nil {
		t.Fatalf("expected no error for valid point, got: %v", err)
	}
	if gotQ != q {
		t.Fatal("expected original query pointer on success")
	}
}

func TestParseIndexesWithGisCrsAndDistanceUnit(t *testing.T) {
	type item struct {
		Pt Point `reindex:"pt,gis,rstar,crs=wgs84,distance_unit=m"`
	}
	joined := map[string][]int{}
	indexes, err := parseIndexes(reflect.TypeOf(item{}), &joined)
	if err != nil {
		t.Fatalf("unexpected parseIndexes error: %v", err)
	}
	if len(indexes) != 1 {
		t.Fatalf("expected 1 index, got %d", len(indexes))
	}
	if indexes[0].IndexType != "gis" {
		t.Fatalf("expected gis index type, got %s", indexes[0].IndexType)
	}
	if indexes[0].Crs != "wgs84" {
		t.Fatalf("expected crs=wgs84, got %s", indexes[0].Crs)
	}
	if indexes[0].DistanceUnit != "m" {
		t.Fatalf("expected distance_unit=m, got %s", indexes[0].DistanceUnit)
	}
}

func TestParseIndexesWithInvalidGisContractOptions(t *testing.T) {
	type badCrs struct {
		Pt Point `reindex:"pt,gis,rstar,crs=epsg4326"`
	}
	type badUnit struct {
		Pt Point `reindex:"pt,gis,rstar,distance_unit=km"`
	}
	type nonGis struct {
		Pt Point `reindex:"pt,rtree,rstar,crs=wgs84"`
	}

	joined := map[string][]int{}
	_, err := parseIndexes(reflect.TypeOf(badCrs{}), &joined)
	if err == nil {
		t.Fatal("expected error for unsupported gis crs")
	}
	_, err = parseIndexes(reflect.TypeOf(badUnit{}), &joined)
	if err == nil {
		t.Fatal("expected error for unsupported gis distance_unit")
	}
	_, err = parseIndexes(reflect.TypeOf(nonGis{}), &joined)
	if err == nil {
		t.Fatal("expected error for crs option on non-gis index")
	}
}

func FuzzGeoDistanceMetersSymmetry(f *testing.F) {
	f.Add(0.0, 0.0, 1.0, 1.0)
	f.Add(-73.9857, 40.7484, 2.2945, 48.8584)
	f.Fuzz(func(t *testing.T, lon1, lat1, lon2, lat2 float64) {
		p1 := Point{lon1, lat1}
		p2 := Point{lon2, lat2}
		if err := ValidateGeoPoint(p1); err != nil {
			t.Skip()
		}
		if err := ValidateGeoPoint(p2); err != nil {
			t.Skip()
		}

		d12 := GeoDistanceMeters(p1, p2)
		d21 := GeoDistanceMeters(p2, p1)
		if math.IsNaN(d12) || math.IsInf(d12, 0) {
			t.Fatalf("invalid distance d12=%v for points %v %v", d12, p1, p2)
		}
		if math.Abs(d12-d21) > 1e-6 {
			t.Fatalf("distance symmetry violation: d12=%v d21=%v", d12, d21)
		}
	})
}

func FuzzGeoDWithinMonotonic(f *testing.F) {
	f.Add(0.0, 0.0, 1.0, 1.0, 1000.0, 2000.0)
	f.Add(-73.9857, 40.7484, -0.1278, 51.5074, 100000.0, 6000000.0)
	f.Fuzz(func(t *testing.T, lon1, lat1, lon2, lat2, r1, r2 float64) {
		p1 := Point{lon1, lat1}
		p2 := Point{lon2, lat2}
		if err := ValidateGeoPoint(p1); err != nil {
			t.Skip()
		}
		if err := ValidateGeoPoint(p2); err != nil {
			t.Skip()
		}
		if r1 < 0 || r2 < 0 {
			t.Skip()
		}
		if r1 > r2 {
			r1, r2 = r2, r1
		}

		withinR1 := GeoDWithin(p1, p2, r1)
		withinR2 := GeoDWithin(p1, p2, r2)
		if withinR1 && !withinR2 {
			t.Fatalf("monotonicity violated: within(%v)=%v within(%v)=%v", r1, withinR1, r2, withinR2)
		}
	})
}

func BenchmarkGeoDistanceMeters(b *testing.B) {
	p1 := Point{-73.9857, 40.7484}
	p2 := Point{-0.1278, 51.5074}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = GeoDistanceMeters(p1, p2)
	}
}

func BenchmarkGeoDWithin(b *testing.B) {
	p1 := Point{-73.9857, 40.7484}
	p2 := Point{-0.1278, 51.5074}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = GeoDWithin(p1, p2, 6_000_000)
	}
}
