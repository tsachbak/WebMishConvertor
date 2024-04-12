using Microsoft.AspNetCore.Mvc;
using WebMishConvertor.Models;

namespace WebMishConvertor.Controllers
{
    [Route("api/[controller]")]
    [ApiController]
    public class ConvertFromITMController : Controller
    {
        [HttpPost]
        [Route("coordinates")]
        public IActionResult ConvertCoordinates([FromBody] ConvertFromItmRequest request)
        {
            int itmEast = request.ITMEast;
            int itmNorth = request.ITMNorth;

            SingleDotPosition dot = new SingleDotPosition(itmEast, itmNorth);

            return Ok(new
            {
                Latitude = dot.Lat,
                Longitude = dot.Longitude,
            });
        }
    }

    public class ConvertFromItmRequest
    {
        public int ITMEast { get; set; }
        public int ITMNorth { get; set; }
    }
}
