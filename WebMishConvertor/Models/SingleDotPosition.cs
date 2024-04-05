namespace WebMishConvertor.Models
{
    public class SingleDotPosition
    {

        private double m_lat;
        private double m_lon;
        private int m_itmEast;
        private int m_itmNorth;

        public SingleDotPosition(double lat, double lon)
        {
            this.m_lat = lat;
            this.m_lon = lon;
            Coordinates.Converters.wgs842itm(m_lat, m_lon, out m_itmNorth, out m_itmEast);
        }

        public SingleDotPosition(int itmNorth, int itmEast)
        {
            this.m_itmEast = itmEast;
            this.m_itmNorth = itmNorth;
            Coordinates.Converters.itm2wgs84(itmNorth, itmEast, out m_lat, out m_lon);
        }

        public double Lat
        {
            get { return m_lat; }
        }

        public double Longitude
        {
            get { return m_lon; }
        }

        public int MitEast
        {
            get { return m_itmEast; }
        }

        public int MitNorth
        {
            get { return m_itmNorth; }
        }
    }
}
