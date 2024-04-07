using Microsoft.AspNetCore.Mvc;
using WebMishConvertor.Models;

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

            SingleDotPosition dot = new SingleDotPosition(latitude, longitude);

            return Ok(new
            {
                ITMEast = dot.MitEast,
                ITMNorth = dot.MitNorth,
            });
        }
    }

    public class CoordinatesRequest
    {
        public double Latitude { get; set;}
        public double Longitude { get; set;}
    }
}
