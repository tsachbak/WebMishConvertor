using Microsoft.AspNetCore.Mvc;

namespace WebMishConvertor.Controllers
{
    [Route("api/[controller]")]
    [ApiController]
    public class ConvertController : Controller
    {
        [HttpPost]
        [Route("coordinates")]
        public IActionResult ConvertCoordinates([FromBody] CoordinatesRequest request)
        {
            double latitude = request.Latitude;
            double longitude = request.Longitude;

            return Ok(new
            {
                ITMEast = latitude,
                ITMNorth = longitude,
            });
        }
    }

    public class CoordinatesRequest
    {
        public double Latitude { get; set;}
        public double Longitude { get; set;}
    }
}
