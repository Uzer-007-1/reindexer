package reindexer

import (
	"fmt"
	"math"
)

const (
	earthRadiusMeters = 6371008.8
	minLongitude      = -180.0
	maxLongitude      = 180.0
	minLatitude       = -90.0
	maxLatitude       = 90.0
)

// ValidateGeoPoint validates a point as WGS84 [lon, lat] pair.
func ValidateGeoPoint(point Point) error {
	if math.IsNaN(point[0]) || math.IsInf(point[0], 0) {
		return fmt.Errorf("invalid longitude value: %v", point[0])
	}
	if math.IsNaN(point[1]) || math.IsInf(point[1], 0) {
		return fmt.Errorf("invalid latitude value: %v", point[1])
	}
	if point[0] < minLongitude || point[0] > maxLongitude {
		return fmt.Errorf("longitude out of range [%v, %v]: %v", minLongitude, maxLongitude, point[0])
	}
	if point[1] < minLatitude || point[1] > maxLatitude {
		return fmt.Errorf("latitude out of range [%v, %v]: %v", minLatitude, maxLatitude, point[1])
	}
	return nil
}

// GeoDistanceMeters returns great-circle distance between 2 WGS84 points in meters.
func GeoDistanceMeters(lhs Point, rhs Point) float64 {
	lat1 := lhs[1] * math.Pi / 180.0
	lat2 := rhs[1] * math.Pi / 180.0
	dLat := lat2 - lat1
	dLon := (rhs[0] - lhs[0]) * math.Pi / 180.0

	sinLat := math.Sin(dLat / 2.0)
	sinLon := math.Sin(dLon / 2.0)
	a := sinLat*sinLat + math.Cos(lat1)*math.Cos(lat2)*sinLon*sinLon
	if a < 0 {
		a = 0
	} else if a > 1 {
		a = 1
	}
	c := 2.0 * math.Atan2(math.Sqrt(a), math.Sqrt(1.0-a))

	return earthRadiusMeters * c
}

// GeoDWithin checks whether two WGS84 points are within specified distance in meters.
func GeoDWithin(lhs Point, rhs Point, distanceMeters float64) bool {
	return GeoDistanceMeters(lhs, rhs) <= distanceMeters
}
